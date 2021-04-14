# Launch a Cromwell server using the Local backend and a local MySQL instance (such as one launched by
# `start_publish_mysql_docker.sh`). Well suited for publishing new versions of Cromwell.
#
# Vault authentication is required to render config files, see
# https://github.com/broadinstitute/dsde-toolbox#authenticating-to-vault.
# Once authenticated, the Cromwell config files can be rendered with:
#
# sbt renderCiResources
#
# Once Cromwell is running, the publish WDL, inputs and options can be submitted using Swagger
# at `http://localhost:8000/swagger/index.html`.

# Fail loudly if things go wrong https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
set -euo pipefail

# https://stackoverflow.com/questions/957928/is-there-a-way-to-get-the-git-root-directory-in-one-command
# shellcheck disable=SC2155
export CROMWELL_ROOT=$(git rev-parse --show-toplevel)

# Set up variables required by the configuration files
export CROMWELL_BUILD_CENTAUR_SLICK_PROFILE=slick.jdbc.MySQLProfile$
export CROMWELL_BUILD_CENTAUR_JDBC_DRIVER=com.mysql.cj.jdbc.Driver
export CROMWELL_BUILD_CENTAUR_JDBC_URL="jdbc:mysql://localhost:3306/cromwell_test?allowPublicKeyRetrieval=true&useSSL=false&rewriteBatchedStatements=true&serverTimezone=UTC&useInformationSchema=true"
export CROMWELL_BUILD_RESOURCES_DIRECTORY="${CROMWELL_ROOT}/target/ci/resources"
export CROMWELL_BUILD_PAPI_JSON_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-service-account.json"
export CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=128000
export CROMWELL_BUILD_CENTAUR_256_BITS_KEY="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA="

# Launch the Cromwell server
cd "${CROMWELL_ROOT}" && sbt -Dconfig.file="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/local_application.conf -Dcall-caching.enabled=true "server/run server"
