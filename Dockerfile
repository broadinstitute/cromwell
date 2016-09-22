# http://github.com/broadinstitute/scala-baseimage
FROM broadinstitute/scala-baseimage

# Cromwell's HTTP Port
EXPOSE 8000

# Install Cromwell
ADD . /cromwell
RUN ["/bin/bash", "-c", "/cromwell/docker/install.sh /cromwell"]

# Add Cromwell as a service (it will start when the container starts)
RUN mkdir /etc/service/cromwell && \
    mkdir /var/log/cromwell && \
    cp /cromwell/docker/run.sh /etc/service/cromwell/run && \
    apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D && \
    echo 'deb https://apt.dockerproject.org/repo ubuntu-trusty main' > /etc/apt/sources.list.d/docker.list && \
    apt-get update && \
    apt-get -yq install docker-engine=1.8.2-0~trusty && \
    apt-get -yq autoremove && \
    apt-get -yq clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /var/tmp/*

# These next 4 commands are for enabling SSH to the container.
# id_rsa.pub is referenced below, but this should be any public key
# that you want to be added to authorized_keys for the root user.
# Copy the public key into this directory because ADD cannot reference
# Files outside of this directory

#EXPOSE 22
#RUN rm -f /etc/service/sshd/down
#ADD id_rsa.pub /tmp/id_rsa.pub
#RUN cat /tmp/id_rsa.pub >> /root/.ssh/authorized_keys
