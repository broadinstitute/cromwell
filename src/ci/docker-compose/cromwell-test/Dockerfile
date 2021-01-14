FROM us.gcr.io/broad-dsp-gcr-public/base/jre:8-debian

RUN apt-get update -q && \
    apt-get upgrade -qq && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /cromwell-test/
COPY docker-setup.sh ./
RUN ./docker-setup.sh
USER hoggett
