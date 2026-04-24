FROM mambaorg/micromamba:latest

COPY environment.yaml .

RUN micromamba create -y -f environment.yaml && \
    micromamba clean --all --yes

COPY . .

ENV ENV_NAME=ChoCallate

WORKDIR /workspace

ENTRYPOINT ["micromamba", "run", "-n", "ChoCallate", "/tmp/main.nf"]
