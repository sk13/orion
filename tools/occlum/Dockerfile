# Orion inside an Occlum SGX enclave: Occlum dev image + Go + conda-forge Python 3.14.
ARG OCCLUM_IMAGE=sk13sk13/occlum:sk13-ubuntu24.04
FROM golang:1.27.3 AS go
FROM ${OCCLUM_IMAGE}
RUN apt-get update && apt-get install -y build-essential libgmp-dev libssl-dev pkg-config && rm -rf /var/lib/apt/lists/*
COPY --from=go /usr/local/go /usr/local/go
ENV GOROOT=/usr/local/go
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/mc.sh && bash /tmp/mc.sh -b -p /opt/miniconda && rm /tmp/mc.sh \
 && /opt/miniconda/bin/conda config --system --remove channels defaults; /opt/miniconda/bin/conda config --system --add channels conda-forge && /opt/miniconda/bin/conda config --system --set channel_priority strict \
 && /opt/miniconda/bin/conda create -y --prefix /opt/python-orion python=3.14 pyyaml numpy scipy matplotlib h5py tqdm certifi pytorch torchvision pip \
 && /opt/miniconda/bin/conda clean -afy && rm -rf /opt/miniconda
