FROM hseeberger/scala-sbt:8u265_1.4.3_2.13.4

RUN apt-get -yq update \
  && apt-get -yq upgrade \
  && apt-get install -yq apt-transport-https ca-certificates curl gnupg-agent \
    jq shellcheck software-properties-common \
  && apt-get -y autoremove \
  && apt-get -yq clean \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /tmp/* \
  && rm -rf /var/tmp/*

# add the docker repo key
RUN curl -fsSL https://download.docker.com/linux/debian/gpg | apt-key add - && \
    curl -fsSL https://apt.releases.hashicorp.com/gpg | apt-key add -
    
# add the docker repo
RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable" && \
    apt-add-repository "deb [arch=amd64] https://apt.releases.hashicorp.com $(lsb_release -cs) main"

RUN apt-get -yq update \
    && apt-get -yq install docker-ce docker-ce-cli containerd.io vault

RUN wget https://github.com/awslabs/git-secrets/archive/1.3.0.zip && \
    unzip 1.3.0.zip && cd git-secrets-1.3.0 && make install && cd .. && rm -rf git-secrets-1.3.0
    
