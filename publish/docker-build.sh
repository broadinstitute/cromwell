#!/usr/bin/env bash

# Builds and pushes the docker image used for publishing cromwell

set -euo pipefail

export DOCKER_CLI_EXPERIMENTAL=enabled

docker buildx rm cromwell-multi-arch-builder || true
docker buildx create --use --name cromwell-multi-arch-builder

build_root="$( dirname "${BASH_SOURCE[0]}" )"
docker buildx build "${build_root}" --platform linux/amd64,linux/arm64 -t broadinstitute/cromwell-publish:latest --push

docker buildx rm cromwell-multi-arch-builder
