FROM broadinstitute/cromwell:develop

RUN git clone https://github.com/vishnubob/wait-for-it.git
RUN mkdir cromwell-working-dir
WORKDIR cromwell-working-dir

COPY ./app-config /app-config

ENTRYPOINT ["/bin/sh", "-c"]
