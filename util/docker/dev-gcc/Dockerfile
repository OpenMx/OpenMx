## Emacs, make this -*- mode: sh; -*-

#-------------- dev-gcc
FROM rocker/r-devel-san as dev-gcc

MAINTAINER "Joshua N. Pritikin" jpritikin@pobox.com

RUN apt-get update && apt-get install -y --no-install-suggests \
			      git \
			      libcurl4-gnutls-dev \
			      libssl-dev \
			      libxml2-dev

WORKDIR /opt/github.org/

RUN git clone -b stable https://github.com/OpenMx/OpenMx.git && \
    cd OpenMx && \
    RD --no-save -f util/update-dependencies.R --args ./DESCRIPTION.in && \
    make cran-install REXEC=RD && \
    rm -rf /opt/github.org/

WORKDIR /root

ENV ASAN_OPTIONS abort_on_error=1

CMD ["R"]

