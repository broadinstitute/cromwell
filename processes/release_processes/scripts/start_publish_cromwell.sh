# Launch a Cromwell server suited for publishing new versions of Cromwell. This uses the Local backend
# and presumes the existence of an appropriately configured MySQL instance. See
# `start_publish_mysql_docker.sh` for launching a containerized version of such a MySQL instance.
#
# Vault authentication is required to render config files, see
# https://github.com/broadinstitute/dsde-toolbox#authenticating-to-vault.
#
# Once Cromwell is running, the publish WDL, inputs and options can be submitted using Swagger
# at `http://localhost:8000/swagger`.
#
# Terminate the Cromwell server and MySQL server (if appropriate to avoid unwanted future call caching)
# when the publish is complete.

# Fail loudly if things go wrong https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
set -euo pipefail

# https://stackoverflow.com/questions/957928/is-there-a-way-to-get-the-git-root-directory-in-one-command
# shellcheck disable=SC2155
export CROMWELL_ROOT=$(git rev-parse --show-toplevel)

cd "${CROMWELL_ROOT}"

sbt renderCiResources

# Set up variables required by the configuration files
export CROMWELL_BUILD_CENTAUR_SLICK_PROFILE=slick.jdbc.MySQLProfile$
export CROMWELL_BUILD_CENTAUR_JDBC_DRIVER=com.mysql.cj.jdbc.Driver
export CROMWELL_BUILD_CENTAUR_JDBC_URL="jdbc:mysql://localhost:3306/cromwell_test?allowPublicKeyRetrieval=true&useSSL=false&rewriteBatchedStatements=true&serverTimezone=UTC&useInformationSchema=true"
export CROMWELL_BUILD_RESOURCES_DIRECTORY="${CROMWELL_ROOT}/target/ci/resources"
export CROMWELL_BUILD_PAPI_JSON_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-service-account.json"
export CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=128000
export CROMWELL_BUILD_CENTAUR_256_BITS_KEY="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA="

# Launch the Cromwell server. Call caching is enabled by default which is helpful to support resuming
# a publish if one of the steps fails transiently.

sbt \
  -Dconfig.file="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/local_application.conf \
  "server/run server"
