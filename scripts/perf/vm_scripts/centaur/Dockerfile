FROM docker-develop:latest

# Run `docker build scripts/docker-develop -t docker-develop` from the cromwell root before building this image
# By pre-compiling centaur we can speed up the startup time of running it in the perf tests
# To rebuild and push this image, run 
# docker build scripts/perf/vm_scripts/centaur -t us.gcr.io/broad-dsde-cromwell-perf/centaur:perf
# docker push us.gcr.io/broad-dsde-cromwell-perf/centaur:perf
#
# Note: because of how memory intensive compiling centaur is, you might need to increase the memory available to docker...
RUN git clone https://github.com/broadinstitute/cromwell.git /code/cromwell

RUN cd /code/cromwell && sbt centaur/it:compile

WORKDIR /code/cromwell
