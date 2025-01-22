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

cromwell::build::slurm::setup_slurm_environment() {
    # Installs the Slurm Workload Manager (WLM) on Ubuntu
    # https://slurm.schedmd.com/
    sudo apt-get update

    sudo apt-get install -y slurm-wlm
    # As an alternative to 'slurm-wlm', you may also try 'slurm-llnl'.
    # If you use 'slurm-llnl', change the below tee path to '/etc/slurm-llnl/slurm.conf'
    # For reasons that are unclear, slurm-llnl doesn't work in Github Actions, and slurm-wlm doesn't work in Travis.
    # We transitioned from slurm-llnl to slurm-wlm when moving from Travis to Github Actions.
    # See the PR associated with ticket WX-888.  

    # Create various directories used by slurm
    sudo mkdir -p /var/run/munge
    sudo mkdir -p /var/spool/slurmd
    sudo chown slurm:slurm /var/spool/slurmd

    # Set up an AppArmor profile for Apptainer to allow non-root users to use it.
    # This became required in Ubuntu 24.04
    sudo tee /etc/apparmor.d/apptainer << 'EOF'
abi <abi/4.0>,
include <tunables/global>

profile apptainer /usr/bin/apptainer{,-suid} flags=(unconfined) {
  userns,
}
EOF

    # Reload to get new profile created above
    sudo systemctl reload apparmor

    # A mash of configure-until-it-runs. Feel free to PR suggestions/fixes.
    # https://slurm.schedmd.com/tutorials.html
    # https://slurm.schedmd.com/configurator.html
    # https://slurm.schedmd.com/slurm.conf.html
    # https://slurm.schedmd.com/quickstart_admin.html
    cat <<SLURM_CONF | sudo tee /etc/slurm/slurm.conf >/dev/null
ClusterName=localhost
ControlMachine=localhost
NodeName=localhost CPUs=4 Sockets=1 CoresPerSocket=2 ThreadsPerCore=2
PartitionName=localpartition Nodes=localhost Default=YES Oversubscribe=Force
ProctrackType=proctrack/pgid
ReturnToService=1
SelectType=select/cons_tres
SelectTypeParameters=CR_CPU
SlurmctldDebug=3
StateSaveLocation=/var/spool/slurmd
SLURM_CONF

    # Start the slurm master
    sudo slurmctld

    # Start the slurm node
    sudo slurmd
}
