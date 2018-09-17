#!/usr/bin/env bash

# Script that runs the perf test using centaur. Invoked in docker.
cd /code/cromwell
git checkout ${CROMWELL_BRANCH}
sbt -Dconfig.file=${PERF_ROOT}/vm_scripts/centaur/centaur.conf "centaur/it:testOnly centaur.ExternalTestCaseSpec"
