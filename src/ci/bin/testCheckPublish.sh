#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::pip_install mkdocs
mkdocs build -s

sbt checkRestApiDocs +package assembly dockerPushCheck +doc

git secrets --scan-history
