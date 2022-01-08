#!/usr/bin/env bash

# Builds and pushes the docker image used for publishing cromwell

set -euo pipefail

build_root="$( dirname "${BASH_SOURCE[0]}" )"
docker build "${build_root}" -t broadinstitute/cromwell-publish
docker push broadinstitute/cromwell-publish
