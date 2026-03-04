FROM ubuntu:24.04 AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install art_modern from GitHub release (.deb for Ubuntu 24.04)
ARG ART_MODERN_VERSION=1.3.4
RUN curl -sL -o /tmp/art-modern.deb \
        "https://github.com/YU-Zhejian/art_modern/releases/download/${ART_MODERN_VERSION}/art-modern_${ART_MODERN_VERSION}%2Bdfsg-1_amd64-ubuntu-2404.deb" \
    && apt-get update \
    && apt-get install -y --no-install-recommends /tmp/art-modern.deb \
    && rm -f /tmp/art-modern.deb \
    && rm -rf /var/lib/apt/lists/*

# Install nocasim
WORKDIR /opt/nocasim
COPY pyproject.toml .
COPY nocasim/ nocasim/
COPY data/ data/

RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir .

ENV PATH="/opt/venv/bin:${PATH}"

ENTRYPOINT ["nocasim"]
CMD ["--help"]
