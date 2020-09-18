FROM openjdk:8

WORKDIR /cromwell-test/
COPY docker-setup.sh ./
RUN ./docker-setup.sh
USER hoggett
ENV PYTHON3_HOME=/usr/bin
