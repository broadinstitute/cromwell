#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::assemble_jars

java -jar $CROMWELL_BUILD_CROMWELL_JAR run ./centaur/src/main/resources/standardTestCases/hello/hello.wdl --inputs ./centaur/src/main/resources/standardTestCases/hello/hello.inputs --metadata-output ./run_mode_metadata.json > console_output.txt

# grep exits 1 if no matches
grep "terminal state: WorkflowSucceededState" console_output.txt
grep "\"wf_hello.hello.salutation\": \"Hello m'Lord!\"" console_output.txt

grep "\"actualWorkflowLanguageVersion\": \"draft-2\"" run_mode_metadata.json
grep "\"actualWorkflowLanguage\": \"WDL\"" run_mode_metadata.json
grep "\"inputs\": \"{\n  \"wf_hello.hello.addressee\\": "m'Lord\"\n}\n\n\"" run_mode_metadata.json
grep "\"wf_hello.hello.salutation\": \"Hello m'Lord!\"" run_mode_metadata.json
grep "\"wf_hello.hello.addressee\": \"Hello m'Lord!\"" run_mode_metadata.json
grep "\"dockerImageUsed\": \"ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950\"" run_mode_metadata.json