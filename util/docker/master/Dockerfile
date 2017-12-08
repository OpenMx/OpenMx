## Emacs, make this -*- mode: sh; -*-

FROM r-base

MAINTAINER "Joshua N. Pritikin" jpritikin@pobox.com

RUN apt-get update && apt-get install -y --no-install-suggests \
			      git \
			      libcurl4-gnutls-dev \
			      libssl-dev \
			      libxml2-dev

WORKDIR /opt/github.org/

RUN git clone -b master https://github.com/OpenMx/OpenMx.git && \
    cd OpenMx && \
    R --no-save -f util/update-dependencies.R --args ./DESCRIPTION.in && \
    make cran-install REXEC=R && \
    rm -rf /opt/github.org/

WORKDIR /root

CMD ["R"]

