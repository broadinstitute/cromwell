#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of common singularity functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::singularity::*
#     Functions for use in other Singularity scripts
#
#   - cromwell::private::singularity::*
#     Functions for use only within this file by cromwell::build::singularity::* functions
#

cromwell::build::singularity::setup_singularity_environment() {
    # Installs one of the singularity forks on Ubuntu
    # https://old.reddit.com/r/HPC/comments/r61bto/singularity_joins_the_linux_foundation_and_is/

    # Using https://sylabs.io/guides/3.9/admin-guide/installation.html#install-from-provided-rpm-deb-packages
    CROMWELL_SINGULARITY_VERSION=3.9.8
    CROMWELL_SINGULARITY_URL="https://github.com/sylabs/singularity/releases/download/v${CROMWELL_SINGULARITY_VERSION}/singularity-ce_${CROMWELL_SINGULARITY_VERSION}-$(lsb_release -cs)_amd64.deb"

    CROMWELL_SINGULARITY_DEB="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/singularity.deb"

    sudo apt-get update

    sudo apt-get install -y squashfs-tools

    curl \
        --location --fail --silent --show-error \
        --output "${CROMWELL_SINGULARITY_DEB}" \
        "${CROMWELL_SINGULARITY_URL}"

    sudo dpkg -i "${CROMWELL_SINGULARITY_DEB}"
}
