## Emacs, make this -*- mode: sh; -*-

MAINTAINER "Joshua N. Pritikin" jpritikin@pobox.com

#-------------- dev-gcc
FROM rocker/r-devel-san as dev-gcc

RUN apt-get update && apt-get install -y --no-install-suggests \
			      git \
			      libcurl4-gnutls-dev \
			      libssl-dev \
			      libxml2-dev

WORKDIR /opt/github.org/

RUN git clone -b master https://github.com/OpenMx/OpenMx.git && \
    cd OpenMx && \
    R --no-save -f util/update-dependencies.R --args ./DESCRIPTION.in && \
    make cran-install && \
    rm -rf /opt/github.org/

WORKDIR /root

ENV ASAN_OPTIONS abort_on_error=1

CMD ["R"]

#-------------- dev-clang
FROM rocker/r-devel-ubsan-clang as dev-clang

RUN apt-get update && apt-get install -y --no-install-suggests \
			      git \
			      libcurl4-gnutls-dev \
			      libssl-dev \
			      libxml2-dev

WORKDIR /opt/github.org/

RUN git clone -b master https://github.com/OpenMx/OpenMx.git && \
    cd OpenMx && \
    R --no-save -f util/update-dependencies.R --args ./DESCRIPTION.in && \
    make cran-install && \
    rm -rf /opt/github.org/

WORKDIR /root

ENV ASAN_OPTIONS abort_on_error=1

CMD ["R"]

#-------------- prod
FROM r-base as prod

RUN apt-get update && apt-get install -y --no-install-suggests \
			      git \
			      libcurl4-gnutls-dev \
			      libssl-dev \
			      libxml2-dev

WORKDIR /opt/github.org/

RUN git clone -b master https://github.com/OpenMx/OpenMx.git && \
    cd OpenMx && \
    R --no-save -f util/update-dependencies.R --args ./DESCRIPTION.in && \
    make cran-install && \
    rm -rf /opt/github.org/

WORKDIR /root

CMD ["R"]

