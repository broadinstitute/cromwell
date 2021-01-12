FROM adoptopenjdk:8-jre-hotspot

RUN apt-get update -q && \
    apt-get upgrade -qq && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /cromwell-test/
COPY docker-setup.sh ./
RUN ./docker-setup.sh
USER hoggett
