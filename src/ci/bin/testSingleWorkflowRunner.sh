#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::start_build_heartbeat

cromwell::build::assemble_jars

java -jar $CROMWELL_BUILD_CROMWELL_JAR run ./centaur/src/main/resources/standardTestCases/hello/hello.wdl --inputs ./centaur/src/main/resources/standardTestCases/hello/hello.inputs --metadata-output ./run_mode_metadata.json | tee console_output.txt

# grep exits 1 if no matches
grep "terminal state: WorkflowSucceededState" console_output.txt
grep "\"wf_hello.hello.salutation\": \"Hello m'Lord!\"" console_output.txt

cat > expected.json <<FIN
{
  "actualWorkflowLanguage": "WDL",
  "actualWorkflowLanguageVersion": "draft-2",
  "dockerImageUsed": "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950",
  "addressee": "m'Lord",
  "salutation": "Hello m'Lord!"
}
FIN

jq '{
  actualWorkflowLanguage,
  actualWorkflowLanguageVersion,
  dockerImageUsed:.calls["wf_hello.hello"][0].dockerImageUsed,
  addressee:.inputs["wf_hello.hello.addressee"],
  salutation:.outputs["wf_hello.hello.salutation"]
}' run_mode_metadata.json > actual.json

cmp <(jq -cS . actual.json) <(jq -cS . expected.json)
