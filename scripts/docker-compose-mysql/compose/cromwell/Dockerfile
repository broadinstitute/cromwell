FROM broadinstitute/cromwell:61-64592db-SNAP

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

RUN git clone https://github.com/vishnubob/wait-for-it.git
RUN mkdir cromwell-working-dir
WORKDIR cromwell-working-dir

COPY ./app-config /app-config

ENTRYPOINT ["/bin/sh", "-c"]
