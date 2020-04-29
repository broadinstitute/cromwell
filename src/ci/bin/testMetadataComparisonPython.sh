#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

build_root="$( dirname "${BASH_SOURCE[0]}" )/../../.."

# Out Python scripts can only be run by Python version >= 3.6, because we use string interpolation which was introduced
# in Python 3.6. But Travis environment is Ubuntu 16.04 Xenial LTS, which supports only Python 3.5.
# Thus we have to run Python tests in Docker container.
docker run -it --rm -v "${build_root}/scripts/metadata_comparison:/metadata_comparison" python:3.6 /bin/bash -c "
pip install --upgrade requests
pip install --upgrade google-api-python-client
pip install --upgrade google-cloud
pip install --upgrade google-cloud-storage
pip install --upgrade pandas
python -m unittest discover -v /metadata_comparison
"
