# Alternatively instead of `FROM linuxbrew/brew` we could run all of the steps used to install brew in docker-setup.sh
# https://github.com/Homebrew/brew/blob/0ff2afdfa8c5943a0e55d9bfe3cdb5d11da8342a/Dockerfile
FROM linuxbrew/brew

WORKDIR /cromwell-publish/
COPY docker-setup.sh git-setup.sh ./
RUN ./docker-setup.sh
