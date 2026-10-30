# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
"""Optional Cython accelerators for ``generate_consensus.py``.

The module is compiled in place (see ``bin/build_consensus_accel.py``).  When it
is not available, ``generate_consensus.py`` falls back to an equivalent
pure-pysam implementation, so the pipeline keeps working from a plain checkout.

``PileupTally``
    Computes FORMAT/DP and the number of reads supporting the reference base for
    the current column of a pysam pileup iterator directly from the htslib
    pileup buffer.  This mirrors ``build_pileup_site``'s overlap handling: one
    observation per query name, highest base quality wins (ties keep the first
    observation), and a deletion is an observation without an allele, so DP
    matches the overlap-corrected mosdepth coverage BED.

``BcfLineWriter``
    Builds records with htslib's ``bcf_update_*`` functions and writes them
    through ``bcf_write``.  Constructing records with
    ``VariantFile.new_record`` plus per-sample assignments costs ~20us/record,
    which dominates the runtime at one record per position; direct htslib
    updates cost a fraction of that.

``ReferenceBlockWriter``
    Drives the reference-only fast path: advances the pileup, tallies the
    reference support without allocating a Python tuple per base, applies the
    absent-caller vote and writes single-base records or contiguous gVCF
    blocks.  This removes the per-base Python loop that otherwise remains even
    when the pileup and record writing are already in C.
"""

from cpython.unicode cimport PyUnicode_AsUTF8AndSize
from libc.stdint cimport int32_t, uint8_t, uint32_t, uint64_t
from libc.stdlib cimport free, malloc
from libc.string cimport memcmp, memset

from pysam.libcalignmentfile cimport IteratorColumn
from pysam.libcbcf cimport VariantHeader
from pysam.libchtslib cimport (
    bam_get_qname,
    bam_get_qual,
    bam_get_seq,
    bam_pileup1_t,
    bam_seqi,
    bcf1_t,
    bcf_clear1,
    bcf_destroy1,
    bcf_hdr_t,
    bcf_hdr_write,
    bcf_init1,
    bcf_update_alleles_str,
    bcf_update_format_int32,
    bcf_update_genotypes,
    bcf_update_info_int32,
    bcf_write,
    hts_close,
    hts_open,
    htsFile,
    seq_nt16_str,
)

DEF MISSING_QUAL = 0xff
# Collapsed-observation state for a deletion: counted in DP, never an allele.
DEF NO_ALLELE = 2
# Table slots above this are refused; the caller falls back to the Python path.
DEF MAX_TABLE_SLOTS = 1 << 21
# The reference fast path grows the table instead of falling back to Python;
# this is only a sanity bound for a pathologically deep column.
DEF FAST_MAX_SLOTS = 1 << 30


cdef inline bytes _alleles_for_base(char base):
    """Pre-joined ``REF,<*>`` allele string for a monomorphic record."""
    if base == 65:  # A
        return b'A,<*>'
    if base == 67:  # C
        return b'C,<*>'
    if base == 71:  # G
        return b'G,<*>'
    if base == 84:  # T
        return b'T,<*>'
    return None


cdef inline uint64_t _hash_qname(const char *name, size_t length) noexcept:
    """FNV-1a hash of a query name (bytes in ``name``)."""
    cdef uint64_t h = 14695981039346656037ULL
    cdef size_t i
    for i in range(length):
        h ^= <uint8_t>name[i]
        h *= 1099511628211ULL
    return h


cdef class PileupTally:
    """Overlap-corrected DP/reference-base counts and full column entries of a
    pysam pileup.

    The object drives the pileup itself (through the iterator's C ``cnext``),
    so no PileupColumn object is created per reference position.  ``sync``
    advances to a target position, ``simple`` tallies the current column and
    ``entries`` materialises it for the rare full-tally path.
    """

    cdef IteratorColumn it
    cdef size_t cap
    cdef uint64_t *hashes
    cdef char **qnames
    cdef uint32_t *name_lens
    cdef uint32_t *best_q
    cdef uint8_t *is_ref
    cdef uint32_t *epoch
    cdef uint32_t cur_epoch
    cdef bint started
    cdef bint done

    def __cinit__(self, object iterator):
        self.it = iterator
        self.cap = 0
        self.hashes = NULL
        self.qnames = NULL
        self.name_lens = NULL
        self.best_q = NULL
        self.is_ref = NULL
        self.epoch = NULL
        self.cur_epoch = 0
        self.started = False
        self.done = False

    def __dealloc__(self):
        free(self.hashes)
        free(self.qnames)
        free(self.name_lens)
        free(self.best_q)
        free(self.is_ref)
        free(self.epoch)

    cdef int _reserve(self, int n) except -1:
        cdef size_t need = <size_t>n * 2 + 8
        cdef size_t new_cap = 64
        if need <= self.cap:
            return 0
        while new_cap < need:
            new_cap <<= 1
        free(self.hashes)
        free(self.qnames)
        free(self.name_lens)
        free(self.best_q)
        free(self.is_ref)
        free(self.epoch)
        self.hashes = <uint64_t *>malloc(new_cap * sizeof(uint64_t))
        self.qnames = <char **>malloc(new_cap * sizeof(char *))
        self.name_lens = <uint32_t *>malloc(new_cap * sizeof(uint32_t))
        self.best_q = <uint32_t *>malloc(new_cap * sizeof(uint32_t))
        self.is_ref = <uint8_t *>malloc(new_cap * sizeof(uint8_t))
        self.epoch = <uint32_t *>malloc(new_cap * sizeof(uint32_t))
        self.cap = 0
        if (self.hashes is NULL or self.qnames is NULL
                or self.name_lens is NULL or self.best_q is NULL
                or self.is_ref is NULL or self.epoch is NULL):
            raise MemoryError()
        memset(self.epoch, 0, new_cap * sizeof(uint32_t))
        self.cap = new_cap
        self.cur_epoch = 0
        return 0

    cpdef bint sync(self, int pos1):
        """Advance the pileup to the first column at or after ``pos1``.

        Returns True if that column is exactly at ``pos1`` (1-based), False if
        it is further ahead or the pileup is exhausted.
        """
        cdef int p
        if self.done:
            return False
        while True:
            if self.started and self.it.pos + 1 >= pos1:
                break
            # bam_mplp_auto returns 1 while columns are available, 0 once the
            # region is exhausted and < 0 on error.
            p = self.it.cnext()
            if p == 0:
                self.done = True
                return False
            if p < 0:
                raise OSError("error while reading the pileup")
            self.started = True
        return self.it.pos + 1 == pos1

    def entries(self):
        """Materialise the current column like pysam's bulk accessors.

        Returns ``(names, sequences, qualities)`` where sequences carry
        samtools-style indel annotations ('A+2TT', 'A-2', '*', '>', '<') and
        qualities use htslib's 255 for missing values, so the shared Python
        full-tally code can consume them unchanged.
        """
        cdef int n = self.it.n_plp
        cdef const bam_pileup1_t *start = self.it.plp
        cdef const bam_pileup1_t *p
        cdef char *qname
        cdef char *buf
        cdef int i, j, l_qname, qpos, indel
        cdef uint32_t qual
        cdef str base, inserted, seq
        names = []
        seqs = []
        quals = []
        for i in range(n):
            p = start + i
            qname = bam_get_qname(p.b)
            l_qname = p.b.core.l_qname
            names.append(qname[:l_qname - 1].decode('utf-8'))
            if p.is_refskip:
                seq = '<' if (p.b.core.flag & 16) else '>'
                qual = 255
            elif p.is_del or p.qpos < 0:
                seq = '*'
                qual = 255
            else:
                qpos = p.qpos
                indel = p.indel
                base = chr(seq_nt16_str[bam_seqi(bam_get_seq(p.b), qpos)])
                if indel > 0:
                    buf = <char *>malloc(indel)
                    if buf is NULL:
                        raise MemoryError()
                    for j in range(indel):
                        buf[j] = seq_nt16_str[bam_seqi(bam_get_seq(p.b),
                                                       qpos + 1 + j)]
                    inserted = buf[:indel].decode('ascii')
                    free(buf)
                    seq = base + '+' + str(indel) + inserted
                elif indel < 0:
                    seq = base + '-' + str(-indel)
                else:
                    seq = base
                qual = (<uint8_t *>bam_get_qual(p.b))[qpos]
            seqs.append(seq)
            quals.append(qual)
        return names, seqs, quals

    def simple(self, str ref_base, int max_slots=MAX_TABLE_SLOTS):
        """Return ``(dp, n_ref, n_obs)`` for the current column, or None if the
        column is too deep for the scratch table (caller falls back to Python).

        ``max_slots`` is the table capacity used for that decision; tests
        lower it to exercise the fallback without a million-read column.

        ``dp`` is the number of collapsed observations and equals ``n_obs``:
        one per distinct query name that is not a reference skip.  This
        matches the overlap-corrected mosdepth coverage BED.  A deletion is an
        observation without an allele, so an overlapping mate's allele always
        wins over it.  ``n_ref`` counts observations whose collapsed allele is
        the plain reference base, so ``dp - n_ref`` is the ``<*>`` depth.
        """
        cdef int dp, n_ref, n_obs, status
        cdef char ref_char = 0
        if ref_base:
            ref_char = <char>ord(ref_base[0])
        status = self._simple_counts(ref_char, max_slots, &dp, &n_ref, &n_obs)
        if status == 1:
            return None
        return (dp, n_ref, n_obs)

    cdef int _simple_counts(self, char ref_char, int max_slots,
                            int* dp_out, int* n_ref_out,
                            int* n_obs_out) except -1:
        """C-level form of ``simple`` used by the reference fast path.

        Returns 0 on success and 1 when the column is deeper than ``max_slots``
        allows.  The count values are written through the output pointers so
        the per-base loop does not allocate a Python tuple.
        """
        cdef int n = self.it.n_plp
        cdef const bam_pileup1_t *start = self.it.plp
        cdef const bam_pileup1_t *p
        cdef uint32_t epoch, qual, prev_q
        cdef uint64_t h
        cdef size_t slot, mask
        cdef char *qname
        cdef uint32_t l_qname
        cdef int qpos, i, n_ref = 0, n_obs = 0, is_match
        cdef int prev_state
        cdef char base

        if n > (max_slots >> 1):
            return 1
        self._reserve(n)
        self.cur_epoch += 1
        if self.cur_epoch == 0:
            memset(self.epoch, 0, self.cap * sizeof(uint32_t))
            self.cur_epoch = 1
        epoch = self.cur_epoch
        mask = self.cap - 1
        for i in range(n):
            p = start + i
            if p.is_refskip:
                continue
            qpos = p.qpos
            if p.is_del or qpos < 0:
                # a deletion is an observation but never an allele: any base
                # or insertion from an overlapping mate replaces it
                is_match = NO_ALLELE
                qual = 0
            else:
                base = seq_nt16_str[bam_seqi(bam_get_seq(p.b), qpos)]
                is_match = 1 if (p.indel == 0 and base == ref_char) else 0
                qual = (<uint8_t *>bam_get_qual(p.b))[qpos]
                if qual == MISSING_QUAL:
                    qual = 0
            qname = bam_get_qname(p.b)
            l_qname = p.b.core.l_qname
            h = _hash_qname(qname, l_qname - 1)
            slot = <size_t>h & mask
            while self.epoch[slot] == epoch:
                if (self.hashes[slot] == h
                        and self.name_lens[slot] == l_qname
                        and memcmp(self.qnames[slot], qname,
                                   l_qname) == 0):
                    prev_state = self.is_ref[slot]
                    prev_q = self.best_q[slot]
                    if prev_state == NO_ALLELE:
                        if is_match != NO_ALLELE:
                            self.is_ref[slot] = is_match
                            self.best_q[slot] = qual
                            n_ref += is_match
                    elif (is_match != NO_ALLELE
                            and qual > prev_q):
                        n_ref -= prev_state
                        self.best_q[slot] = qual
                        self.is_ref[slot] = is_match
                        n_ref += is_match
                    break
                slot = (slot + 1) & mask
            else:
                self.epoch[slot] = epoch
                self.hashes[slot] = h
                self.qnames[slot] = qname
                self.name_lens[slot] = l_qname
                self.best_q[slot] = qual
                self.is_ref[slot] = is_match
                if is_match != NO_ALLELE:
                    n_ref += is_match
                n_obs += 1
        dp_out[0] = n_obs
        n_ref_out[0] = n_ref
        n_obs_out[0] = n_obs
        return 0


cdef class BcfLineWriter:
    """Write consensus records to a BCF file through htslib.

    Records are built with htslib's ``bcf_update_*`` functions instead of
    pysam ``VariantRecord`` objects (``new_record`` + per-sample assignments
    costs ~20us per record, which dominates runtime at one record per
    position).
    """

    cdef htsFile *fp
    cdef bcf_hdr_t *hdr
    cdef bcf1_t *rec
    cdef int32_t gt_buf[16]
    cdef int32_t ad_buf[64]
    cdef bint closed

    def __cinit__(self, VariantHeader header, str path, str mode="wb"):
        cdef bytes bpath
        cdef bytes bmode
        self.fp = NULL
        self.rec = NULL
        self.closed = True
        self.hdr = header.ptr
        if self.hdr is NULL:
            raise ValueError("header is not initialized")
        header._hdr_sync()
        bpath = path.encode("utf-8")
        bmode = mode.encode("ascii")
        self.fp = hts_open(<const char *>bpath, <const char *>bmode)
        if self.fp is NULL:
            raise OSError("could not open '%s' for writing" % path)
        if bcf_hdr_write(self.fp, self.hdr) < 0:
            hts_close(self.fp)
            self.fp = NULL
            raise OSError("could not write header to '%s'" % path)
        self.rec = bcf_init1()
        if self.rec is NULL:
            raise MemoryError("could not allocate BCF record")
        self.closed = False

    cdef int _write(self, int rid, int pos0, const char *alleles,
                    const int32_t *gt, int n_gt, const int32_t *ad,
                    int n_ad, int dp, int end_pos=-1) except -1:
        cdef int32_t dp_value = dp
        cdef int32_t end_value
        cdef bcf1_t *rec = self.rec
        bcf_clear1(rec)
        rec.rid = rid
        rec.pos = pos0
        if bcf_update_alleles_str(self.hdr, rec, alleles) < 0:
            raise MemoryError("could not set alleles")
        if end_pos >= 0:
            # symbolic alleles carry INFO/END (1-based, inclusive)
            end_value = end_pos
            if bcf_update_info_int32(self.hdr, rec, b"END", &end_value, 1) < 0:
                raise MemoryError("could not set END")
            # rlen drives CSI overlap queries; without it an indexed multi-base
            # block is invisible to fetches that start inside the block.
            rec.rlen = end_pos - pos0
        if n_gt > 0 and bcf_update_genotypes(self.hdr, rec, gt, n_gt) < 0:
            raise MemoryError("could not set GT")
        if n_ad > 0 and bcf_update_format_int32(self.hdr, rec, b"AD", ad,
                                                n_ad) < 0:
            raise MemoryError("could not set AD")
        if bcf_update_format_int32(self.hdr, rec, b"DP", &dp_value, 1) < 0:
            raise MemoryError("could not set DP")
        if bcf_write(self.fp, self.hdr, rec) < 0:
            raise OSError("could not write record")
        return 0

    cdef int _write_monomorphic(self, int rid, int pos0, bytes alleles,
                                int ad_ref, int ad_other, int dp,
                                int end_pos) except -1:
        """C-level writer used by both the public method and the fast path."""
        cdef int32_t ad[2]
        if end_pos < 0:
            end_pos = pos0 + 1
        self.gt_buf[0] = 2  # bcf_gt_unphased(0)
        self.gt_buf[1] = 2
        ad[0] = ad_ref
        ad[1] = ad_other
        return self._write(rid, pos0, <const char *>alleles, self.gt_buf, 2,
                           ad, 2, dp, end_pos)

    def write_monomorphic(self, int rid, int pos0, bytes alleles,
                          int ad_ref, int ad_other, int dp, int end_pos=-1):
        """Write a gVCF-style reference record or block ``alleles``
        (precomputed, e.g. b'A,<*>'), GT 0/0, AD (ad_ref, ad_other) and DP.

        ``ad_other`` is the depth of reads that do not support the reference
        base (the ``<*>`` slot).  ``end_pos`` is the 1-based inclusive end; a
        negative value means the single base at ``pos0``.
        """
        return self._write_monomorphic(rid, pos0, alleles, ad_ref, ad_other,
                                       dp, end_pos)

    cdef int _write_reference_block(self, int rid, int pos0, bytes alleles,
                                    int dp, int end_pos) except -1:
        """C-level writer for a multi-base gVCF reference block."""
        self.gt_buf[0] = 2  # bcf_gt_unphased(0)
        self.gt_buf[1] = 2
        return self._write(rid, pos0, <const char *>alleles, self.gt_buf, 2,
                           NULL, 0, dp, end_pos)

    def write_reference_block(self, int rid, int pos0, bytes alleles, int dp,
                              int end_pos):
        """Write a multi-base gVCF reference block without FORMAT/AD.

        A single AD value cannot describe every position of the block, and
        ``bcftools merge --gvcf`` reports missing depth when the block is used
        to fill a hom-reference sample at another sample's variant.
        """
        return self._write_reference_block(rid, pos0, alleles, dp, end_pos)

    def write_record(self, int rid, int pos0, object alleles, object gt,
                     object ad, int dp):
        """Write an arbitrary consensus record.

        ``alleles`` is the REF allele followed by ALTs (sequence of str or
        pre-joined bytes), ``gt`` a sequence of allele indices (None for a
        missing allele) and ``ad`` a sequence of per-allele depths.
        """
        cdef int32_t *gt_ptr = self.gt_buf
        cdef int32_t *ad_ptr = self.ad_buf
        cdef int32_t *heap = NULL
        cdef bytes alleles_bytes
        cdef object value
        cdef int i, n_gt, n_ad
        if isinstance(alleles, bytes):
            alleles_bytes = alleles
        else:
            alleles_bytes = ",".join(alleles).encode("ascii")
        gt_ = list(gt)
        ad_ = list(ad)
        n_gt = len(gt_)
        n_ad = len(ad_)
        if n_gt > 16 or n_ad > 64:
            heap = <int32_t *>malloc((n_gt + n_ad) * sizeof(int32_t))
            if heap is NULL:
                raise MemoryError()
            gt_ptr = heap
            ad_ptr = heap + n_gt
        try:
            for i in range(n_gt):
                value = gt_[i]
                gt_ptr[i] = 0 if value is None else ((<int>value + 1) << 1)
            for i in range(n_ad):
                ad_ptr[i] = <int32_t>ad_[i]
            return self._write(rid, pos0, <const char *>alleles_bytes, gt_ptr,
                               n_gt, ad_ptr, n_ad, dp)
        finally:
            free(heap)

    def close(self):
        cdef htsFile *fp = self.fp
        self.fp = NULL
        if self.rec is not NULL:
            bcf_destroy1(self.rec)
            self.rec = NULL
        self.closed = True
        # hts_close flushes the final BGZF block; a failure (e.g. disk full)
        # would otherwise leave a truncated BCF and a successful exit status.
        if fp is not NULL and hts_close(fp) < 0:
            raise OSError("error closing BCF output")


cdef class ReferenceBlockWriter:
    """Write runs of hom-reference positions without entering Python per base.

    The Python driver still decides which positions have caller records and
    where a run starts and ends.  This class walks the BED run, advances the
    pileup, tallies reference support, applies the absent-caller vote and
    emits either single-base records or contiguous gVCF blocks.  It keeps a
    pending block across calls so adjacent BED intervals are not split.
    """

    cdef BcfLineWriter writer
    cdef int threshold
    cdef bint absent_vote_ref
    cdef int min_coverage
    cdef bint compress_blocks
    cdef int n_files
    cdef bint pending
    cdef int pending_rid
    cdef int pending_start0
    cdef int pending_end0
    cdef int pending_dp
    cdef bint pending_ge
    cdef char pending_ref
    cdef int pending_ad_ref
    cdef int pending_ad_other

    def __cinit__(self, BcfLineWriter writer, int threshold,
                  bint absent_vote_ref, int min_coverage,
                  bint compress_blocks, int n_files):
        if writer is None:
            raise ValueError("writer is required")
        self.writer = writer
        self.threshold = threshold
        self.absent_vote_ref = absent_vote_ref
        self.min_coverage = min_coverage
        self.compress_blocks = compress_blocks
        self.n_files = n_files
        self.pending = False

    cdef int _flush(self) except -1:
        """Emit the pending block, if any."""
        cdef bytes alleles
        if not self.pending:
            return 0
        alleles = _alleles_for_base(self.pending_ref)
        if self.pending_end0 != self.pending_start0 + 1:
            self.writer._write_reference_block(
                self.pending_rid, self.pending_start0, alleles,
                self.pending_dp, self.pending_end0)
        else:
            self.writer._write_monomorphic(
                self.pending_rid, self.pending_start0, alleles,
                self.pending_ad_ref, self.pending_ad_other,
                self.pending_dp, self.pending_end0)
        self.pending = False
        return 0

    def flush(self):
        """Write a pending reference block, if one is open."""
        return self._flush()

    def write_interval(self, PileupTally tally, int rid, str ref_seq,
                       int ref_offset, int start1, int end1,
                       int deletion_end):
        """Write the reference-only positions ``start1..end1`` (1-based).

        ``ref_seq`` is the reference sequence of the enclosing BED interval;
        ``ref_offset`` is the 0-based offset of ``start1`` within it.  No
        caller record may fall inside the interval.  ``deletion_end`` is the
        1-based inclusive end of the furthest called deletion, whose REF span
        suppresses reference records.
        """
        cdef Py_ssize_t seq_len
        cdef const char* seq
        cdef long long pos
        cdef int dp, n_ref, n_obs, ref_votes, status
        cdef char ref_base
        cdef bint has_column, ge
        if start1 > end1:
            return
        seq = PyUnicode_AsUTF8AndSize(ref_seq, &seq_len)
        if seq is NULL:
            raise ValueError("reference sequence is not valid UTF-8")
        if (ref_offset < 0 or start1 < 1 or end1 < start1
                or ref_offset + (end1 - start1 + 1) > seq_len):
            raise ValueError("reference sequence does not cover the interval")
        # The loop counter is 64-bit so a BED interval ending at 2**31 - 1
        # does not overflow when the stop value is incremented.
        for pos in range(<long long>start1, <long long>end1 + 1):
            if pos <= deletion_end:
                # The called deletion's REF span already represents this base.
                self._flush()
                continue
            ref_base = seq[ref_offset + pos - start1]
            if not (ref_base == 65 or ref_base == 67
                    or ref_base == 71 or ref_base == 84):
                # IUPAC ambiguity codes (N included) cannot be written.
                self._flush()
                continue
            has_column = tally.sync(<int>pos)
            if has_column:
                status = tally._simple_counts(
                    ref_base, FAST_MAX_SLOTS, &dp, &n_ref, &n_obs)
                if status != 0:
                    raise MemoryError("pileup column is too deep")
            else:
                dp = 0
                n_ref = 0
                n_obs = 0
            if self.absent_vote_ref:
                ref_votes = self.n_files
            elif has_column and dp >= self.min_coverage:
                ref_votes = self.n_files
            else:
                ref_votes = 0
            if ref_votes < self.threshold:
                self._flush()
                continue
            ge = dp >= self.min_coverage
            if not self.compress_blocks:
                self.writer._write_monomorphic(
                    rid, <int>pos - 1, _alleles_for_base(ref_base), n_ref,
                    n_obs - n_ref, dp, <int>pos)
                continue
            if (self.pending and self.pending_rid == rid
                    and <int>pos - 1 == self.pending_end0
                    and ge == self.pending_ge):
                self.pending_end0 = <int>pos
                if dp < self.pending_dp:
                    self.pending_dp = dp
            else:
                self._flush()
                self.pending = True
                self.pending_rid = rid
                self.pending_start0 = <int>pos - 1
                self.pending_end0 = <int>pos
                self.pending_dp = dp
                self.pending_ge = ge
                self.pending_ref = ref_base
                self.pending_ad_ref = n_ref
                self.pending_ad_other = n_obs - n_ref
