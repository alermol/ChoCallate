FROM mambaorg/micromamba:latest

COPY environment.yaml .

RUN micromamba create -y -f environment.yaml && \
    micromamba clean --all --yes

COPY . .

# Build the optional Cython accelerator for generate_consensus.py in place.
# generate_consensus.py falls back to an equivalent pure-pysam implementation if
# the module is missing, so this is an optimisation, not a requirement.
RUN micromamba run -n ChoCallate python bin/build_consensus_accel.py

ENV ENV_NAME=ChoCallate

WORKDIR /workspace

ENTRYPOINT ["micromamba", "run", "-n", "ChoCallate", "/tmp/main.nf"]
