FROM broadinstitute/cromwell:25

RUN apt-get update
RUN apt-get --assume-yes dist-upgrade
RUN apt-get --assume-yes install mysql-client

COPY ./app-setup /app-setup
COPY ./app-config /app-config

WORKDIR /app-setup

RUN chmod u+x robustStart.sh

ENTRYPOINT ["/app-setup/robustStart.sh", "java -Dconfig.file=/app-config/application.conf -jar /app/cromwell.jar server"]
