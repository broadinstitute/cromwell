FROM broadinstitute/cromwell:87-e6e1ebd-SNAP

RUN git clone https://github.com/vishnubob/wait-for-it.git
RUN mkdir cromwell-working-dir
WORKDIR cromwell-working-dir

COPY ./app-config /app-config

ENTRYPOINT ["/bin/sh", "-c"]
