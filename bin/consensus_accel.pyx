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
    observation).

``BcfLineWriter``
    Writes records from VCF text lines through htslib's ``vcf_parse`` +
    ``bcf_write``.  Constructing records with ``VariantFile.new_record`` plus
    per-sample assignments costs ~20us/record, which dominates the runtime at
    one record per position; text lines cost a fraction of that.
"""

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
# Table slots above this are refused; the caller falls back to the Python path.
DEF MAX_TABLE_SLOTS = 1 << 21


cdef inline uint64_t _hash_qname(const char *name, size_t length) noexcept:
    """FNV-1a hash of a query name (bytes in ``name``)."""
    cdef uint64_t h = 14695981039346656037ULL
    cdef size_t i
    for i in range(length):
        h ^= <uint8_t>name[i]
        h *= 1099511628211ULL
    return h


cdef class PileupTally:
    """DP/reference-base counts and full column entries of a pysam pileup.

    The object drives the pileup itself (through the iterator's C ``cnext``),
    so no PileupColumn object is created per reference position.  ``sync``
    advances to a target position, ``simple`` tallies the current column and
    ``entries`` materialises it for the rare full-tally path.
    """

    cdef IteratorColumn it
    cdef size_t cap
    cdef uint64_t *hashes
    cdef char **qnames
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
        self.best_q = NULL
        self.is_ref = NULL
        self.epoch = NULL
        self.cur_epoch = 0
        self.started = False
        self.done = False

    def __dealloc__(self):
        free(self.hashes)
        free(self.qnames)
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
        free(self.best_q)
        free(self.is_ref)
        free(self.epoch)
        self.hashes = <uint64_t *>malloc(new_cap * sizeof(uint64_t))
        self.qnames = <char **>malloc(new_cap * sizeof(char *))
        self.best_q = <uint32_t *>malloc(new_cap * sizeof(uint32_t))
        self.is_ref = <uint8_t *>malloc(new_cap * sizeof(uint8_t))
        self.epoch = <uint32_t *>malloc(new_cap * sizeof(uint32_t))
        self.cap = 0
        if (self.hashes is NULL or self.qnames is NULL or self.best_q is NULL
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

        Returns ``(names, sequences, qualities, dp)`` where sequences carry
        samtools-style indel annotations ('A+2TT', 'A-2', '*', '>', '<') and
        qualities use htslib's 255 for missing values, so the shared Python
        full-tally code can consume them unchanged.
        """
        cdef int n = self.it.n_plp
        cdef const bam_pileup1_t *start = self.it.plp
        cdef const bam_pileup1_t *p
        cdef char *qname
        cdef char *buf
        cdef int i, j, l_qname, qpos, dp = 0, indel
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
                dp += 1
                seq = '*'
                qual = 255
            else:
                dp += 1
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
        return names, seqs, quals, dp

    def simple(self, str ref_base):
        """Return ``(dp, n_ref, n_obs)`` for the current column, or None if the
        column is too deep for the scratch table (caller falls back to Python).

        ``dp`` counts every covering read except reference skips (matching
        ``samtools depth -J``).  ``n_ref`` counts reads whose collapsed allele
        is the plain reference base.  ``n_obs`` counts distinct query names with
        an allele at the column (base, insertion or deletion); ``n_obs - n_ref``
        is the ``<*>`` (other-allele) depth.
        """
        cdef int n = self.it.n_plp
        cdef const bam_pileup1_t *start = self.it.plp
        cdef const bam_pileup1_t *p
        cdef char ref_char
        cdef uint32_t epoch, qual, prev_q
        cdef uint64_t h
        cdef size_t slot, mask
        cdef char *qname
        cdef int l_qname, qpos, i, dp = 0, n_ref = 0, n_obs = 0, is_match
        cdef char base

        if n > (MAX_TABLE_SLOTS >> 1):
            return None
        if ref_base:
            ref_char = <char>ord(ref_base[0])
        else:
            ref_char = 0
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
            dp += 1
            qpos = p.qpos
            if p.is_del or qpos < 0:
                # a deletion at this position is an observation that does not
                # support the reference base
                is_match = 0
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
                        and memcmp(self.qnames[slot], qname, l_qname) == 0):
                    prev_q = self.best_q[slot]
                    if qual > prev_q:
                        n_ref -= self.is_ref[slot]
                        self.best_q[slot] = qual
                        self.is_ref[slot] = is_match
                        n_ref += is_match
                    break
                slot = (slot + 1) & mask
            else:
                self.epoch[slot] = epoch
                self.hashes[slot] = h
                self.qnames[slot] = qname
                self.best_q[slot] = qual
                self.is_ref[slot] = is_match
                n_ref += is_match
                n_obs += 1
        return (dp, n_ref, n_obs)


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

    def __cinit__(self, VariantHeader header, str path):
        cdef bytes bpath
        self.fp = NULL
        self.rec = NULL
        self.closed = True
        self.hdr = header.ptr
        if self.hdr is NULL:
            raise ValueError("header is not initialized")
        header._hdr_sync()
        bpath = path.encode("utf-8")
        self.fp = hts_open(<const char *>bpath, b"wb")
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

    def write_monomorphic(self, int rid, int pos0, bytes alleles,
                          int ad_ref, int ad_other, int dp, int end_pos=-1):
        """Write a gVCF-style reference record or block ``alleles``
        (precomputed, e.g. b'A,<*>'), GT 0/0, AD (ad_ref, ad_other) and DP.

        ``ad_other`` is the depth of reads that do not support the reference
        base (the ``<*>`` slot).  ``end_pos`` is the 1-based inclusive end; a
        negative value means the single base at ``pos0``.
        """
        cdef int32_t ad[2]
        if end_pos < 0:
            end_pos = pos0 + 1
        self.gt_buf[0] = 2  # bcf_gt_unphased(0)
        self.gt_buf[1] = 2
        ad[0] = ad_ref
        ad[1] = ad_other
        return self._write(rid, pos0, <const char *>alleles, self.gt_buf, 2,
                           ad, 2, dp, end_pos)

    def write_reference_block(self, int rid, int pos0, bytes alleles, int dp,
                              int end_pos):
        """Write a multi-base gVCF reference block without FORMAT/AD.

        A single AD value cannot describe every position of the block, and
        ``bcftools merge --gvcf`` reports missing depth when the block is used
        to fill a hom-reference sample at another sample's variant.
        """
        self.gt_buf[0] = 2  # bcf_gt_unphased(0)
        self.gt_buf[1] = 2
        return self._write(rid, pos0, <const char *>alleles, self.gt_buf, 2,
                           NULL, 0, dp, end_pos)

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
        if self.fp is not NULL:
            hts_close(self.fp)
            self.fp = NULL
        if self.rec is not NULL:
            bcf_destroy1(self.rec)
            self.rec = NULL
        self.closed = True
