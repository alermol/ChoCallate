FROM mambaorg/micromamba:latest

COPY environment.yaml .

RUN micromamba create -y -f environment.yaml && \
    micromamba clean --all --yes

# Install the pipeline outside /tmp. Apptainer and Singularity bind the host
# /tmp into the container by default, which would hide a pipeline installed
# there.
#
# The COPY destination is created with the ownership given by --chown, so the
# image user can create NXF_HOME below it. The image runs as MAMBA_USER and
# the accelerator build writes into bin/, so the copied project must be owned
# by that user as well.
COPY --chown=${MAMBA_USER}:${MAMBA_USER} . /opt/chocallate/
WORKDIR /opt/chocallate

# Build the optional Cython accelerator for generate_consensus.py in place.
# generate_consensus.py falls back to an equivalent pure-pysam implementation
# if the module is missing, so this is an optimisation, not a requirement.
RUN micromamba run -n ChoCallate python bin/build_consensus_accel.py

# Pin the plugins in the image so a run does not need network access.
ENV NXF_HOME=/opt/chocallate/.nextflow
RUN micromamba run -n ChoCallate nextflow plugin install \
        nf-schema@3.0.0,nf-boost@0.7.0 && \
    chmod -R a+rwX "${NXF_HOME}"

# Both micromamba run and Nextflow need a writable HOME. With docker run -u
# $(id -u) the UID has no passwd entry, so ~ is not writable. Apptainer and
# Singularity bind the host /tmp, which is why HOME is not set to /tmp either.
ENV HOME=/opt/chocallate/.home
RUN mkdir -p "${HOME}" && chmod 0777 "${HOME}"

ENV ENV_NAME=ChoCallate

WORKDIR /workspace

ENTRYPOINT ["micromamba", "run", "-n", "ChoCallate", \
            "nextflow", "run", "/opt/chocallate/main.nf"]
