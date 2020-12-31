#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

export CROMWELL_BUILD_REQUIRES_SECURE=true

# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

# Our Python scripts can only be run by Python version >= 3.6, because we use string interpolation which was introduced
# in Python 3.6. But Travis environment is Ubuntu 16.04 Xenial LTS, which officially supports only Python <= 3.5.
# Thus we have to run Python tests in Docker container.
GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON="cromwell-centaur-service-account.json"
export GOOGLE_APPLICATION_CREDENTIALS="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON}"
export DIGESTER_TEST_GCS=true

docker run --rm \
           -e GOOGLE_APPLICATION_CREDENTIALS \
           -e DIGESTER_TEST_GCS \
           -v "${CROMWELL_BUILD_RESOURCES_DIRECTORY}:${CROMWELL_BUILD_RESOURCES_DIRECTORY}" \
           -v "${CROMWELL_BUILD_ROOT_DIRECTORY}/scripts/metadata_comparison:/metadata_comparison" \
           python:3 /bin/bash -c "

pip install --upgrade requests
pip install --upgrade google-api-python-client
pip install --upgrade google-cloud
pip install --upgrade google-cloud-storage
pip install --upgrade gitpython
pip install --upgrade python-dateutil
python -m unittest discover -v /metadata_comparison
"
