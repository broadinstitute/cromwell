set -euo pipefail

# Launch a Cromwell server using the Local backend and configured to connect to a locally running MySQL
# instance such as one launched by `start_publish_mysql_docker.sh`. Well suited for publishing new
# versions of Cromwell.
#
# To run this script the current working directory should be the root of the Cromwell git working tree.
# Vault authentication is required to render config files, see
# https://github.com/broadinstitute/dsde-toolbox#authenticating-to-vault.
# Once authenticated, the Cromwell config files can be rendered with:
#
# sbt renderCiResources
#
# Once Cromwell is running, one can submit the publish WDL, inputs and options using Swagger
# at `http://localhost:8000/swagger/index.html`.

# Set up a number of variables required by the configuration files:
export CROMWELL_BUILD_CENTAUR_SLICK_PROFILE=slick.jdbc.MySQLProfile$
export CROMWELL_BUILD_CENTAUR_JDBC_DRIVER=com.mysql.cj.jdbc.Driver
export CROMWELL_BUILD_CENTAUR_JDBC_URL="jdbc:mysql://localhost:3306/cromwell_test?allowPublicKeyRetrieval=true&useSSL=false&rewriteBatchedStatements=true&serverTimezone=UTC&useInformationSchema=true"
export CROMWELL_BUILD_RESOURCES_DIRECTORY="${PWD}/target/ci/resources"
export CROMWELL_BUILD_PAPI_JSON_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-service-account.json"
export CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=128000
export CROMWELL_BUILD_CENTAUR_256_BITS_KEY="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA="

sbt -Dconfig.file=${CROMWELL_BUILD_RESOURCES_DIRECTORY}/local_application.conf "server/run server"