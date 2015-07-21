# http://github.com/broadinstitute/scala-baseimage
FROM broadinstitute/scala-baseimage

# Cromwell's HTTP Port
EXPOSE 8000

# Install Cromwell
ADD . /cromwell
RUN ["/bin/bash", "-c", "/cromwell/docker/install.sh /cromwell"]

# Add Cromwell as a service (it will start when the container starts)
RUN mkdir /etc/service/cromwell && \
    cp /cromwell/docker/run.sh /etc/service/cromwell/run

# These next 4 commands are for enabling SSH to the container.
# id_rsa.pub is referenced below, but this should be any public key
# that you want to be added to authorized_keys for the root user.
# Copy the public key into this directory because ADD cannot reference
# Files outside of this directory

#EXPOSE 22
#RUN rm -f /etc/service/sshd/down
#ADD id_rsa.pub /tmp/id_rsa.pub
#RUN cat /tmp/id_rsa.pub >> /root/.ssh/authorized_keys
