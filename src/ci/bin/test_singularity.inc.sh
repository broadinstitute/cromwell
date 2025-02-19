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

    # Using https://apptainer.org/docs/admin/1.3/installation.html#install-ubuntu-packages
    sudo apt install -y software-properties-common
    sudo add-apt-repository -y ppa:apptainer/ppa
    sudo apt update
    sudo apt install -y apptainer
}
