#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

export AWS_SHARED_CREDENTIALS_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_credentials
export AWS_CONFIG_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_config

# "space" and "scatter" tests are disabled because hey intermittently fail on AWS
# https://broadworkbench.atlassian.net/browse/BA-6152
# NOTE! For some reason, exclusions must all be lower case even if the test name is a mixture of upper and lower

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!! https://broadworkbench.atlassian.net/browse/BT-57 !!!!!!!!!!
# !!!!!!!!!! DISABLING ALL TESTS BUT HELLO WORLD (Travis) AND  !!!!!!!!!!
# !!!!!!!!!!   MUTECT2 (Jenkins) UNTIL AWS CAN AUTHENTICATE    !!!!!!!!!!
# !!!!!!!!!!                WITH DOCKER REPOS                  !!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cromwell::build::run_centaur \
     -i hello \
     -i mutect2.aws

#    -p 100 \
#    -e smartseq2singlesample \
#    -e arrays \
#    -e haplotypecaller \
#    -e jointdiscovery \
#    -e mutect2 \
#    -e singlesample \
#    -e singlesample_production \
#    -e cnv_somatic_pair \
#    -e cnv_somatic_panel \
#    -e localdockertest \
#    -e non_root_default_user \
#    -e non_root_specified_user \
#    -e abort.scheduled_abort \
#    -e relative_output_paths \
#    -e relative_output_paths_colliding \
#    -e standard_output_paths_colliding_prevented \
#    -e space \
#    -e scatter \
#    -e long_cmd \
#    -e runtwiceexpectingcallcaching \
#    -e cwl_cache_within_workflow \
#    -e cachewithinwf

cromwell::build::generate_code_coverage
