#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

build_root="$( dirname "${BASH_SOURCE[0]}" )/../../.."
docker run -it --rm -v "${build_root}/scripts/metadata_comparison:/metadata_comparison" python:3.6 /bin/bash -c "
pip install --upgrade requests
pip install --upgrade google-api-python-client
pip install --upgrade google-cloud
pip install --upgrade google-cloud-storage
pip install --upgrade pandas
python -m unittest discover -v /metadata_comparison
"

cromwell::build::publish_artifacts
