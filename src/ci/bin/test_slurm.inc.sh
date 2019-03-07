#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of common SLURM functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::slurm::*
#     Functions for use in other SLURM scripts
#
#   - cromwell::private::slurm::*
#     Functions for use only within this file by cromwell::build::slurm::* functions
#

cromwell::build::slurm::slurm_install() {
    # Installs the Slurm Workload Manager (WLM) on Ubuntu
    # https://slurm.schedmd.com/
    sudo apt-get update

    # Try the Lawrence Livermore National Laboratory (LLNL) version first
    sudo apt-get install -y slurm-llnl || apt-get install -y slurm-wlm

    # Create various directories used by slurm
    sudo mkdir -p /etc/slurm-llnl
    sudo mkdir -p /var/run/slurm-llnl
    sudo mkdir -p /var/run/munge
    sudo mkdir -p /var/spool/slurmd

    # A mash of configure-until-it-runs. Feel free to PR suggestions/fixes
    # https://slurm.schedmd.com/tutorials.html
    # https://slurm.schedmd.com/configurator.html
    # https://slurm.schedmd.com/slurm.conf.html
    # https://slurm.schedmd.com/quickstart_admin.html
    cat <<SLURM_CONF | sudo tee /etc/slurm-llnl/slurm.conf >/dev/null
ControlMachine=localhost
NodeName=localhost
PartitionName=localpartition Nodes=localhost Default=YES
ProctrackType=proctrack/pgid
ReturnToService=1
SelectType=select/cons_res
SelectTypeParameters=CR_CPU
SlurmctldDebug=3
SLURM_CONF

    # Munge will run as the munge user, so change the owner of the directory
    sudo chown munge:munge /var/run/munge

    # Create a munge key, required for slurm-llnl
    dd if=/dev/random bs=1 count=1024 | sudo tee /etc/munge/munge.key >/dev/null
    sudo chown munge:munge /etc/munge/munge.key
    sudo chmod 400 /etc/munge/munge.key

    # Start munge as the munge user
    sudo -u munge munged

    # Start the slurm master
    sudo slurmctld

    # Start the slurm node
    sudo slurmd
}
