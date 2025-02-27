# This script will launch a Postgres 11.3 container with parameters matching the expectations of Cromwell's CI configuration
# You can run cromwell locally using this container as its DB.
# Your intelliJ repo template or cromwell config should provide the following environment variables:
# CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="org.postgresql.Driver"
# CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="slick.jdbc.PostgresProfile$"
# CROMWELL_BUILD_CENTAUR_JDBC_URL="jdbc:postgresql://localhost:5432/cromwell_test?reWriteBatchedInserts=true"

set -euo pipefail

CURRENT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

docker run \
    --name=cromwell_postgres_test_db \
    --env POSTGRES_USER=cromwell \
    --env POSTGRES_PASSWORD=test \
    --env POSTGRES_DB=cromwell_test \
    --volume ${CURRENT_DIR}/postgresql-initdb.d:/docker-entrypoint-initdb.d \
    --publish 5432:5432 \
    --rm postgres:14
