FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    gromacs \
    python3 \
    python3-pip \
    curl \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages \
    MDAnalysis \
    pandas \
    numpy

WORKDIR /workspace
COPY scripts/ ./scripts/
COPY config/  ./config/
RUN chmod +x scripts/run_all.sh
