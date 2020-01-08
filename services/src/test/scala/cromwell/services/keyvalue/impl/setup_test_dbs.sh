#!/bin/bash

# setup_test_dbs.sh:
# The tests in KeyValueDatabaseSpec assume a number of test databases will be running locally. The
# tests print out a message on how to set these up: This script just runs all of them in sequence
# with some retries in case the docker images aren't available when the DB reset is run.
#
# WARNING: Run this from the cromwell checkout root, or the directory links will fail.
# WARNING2: If the setup text in the tests change, this script will be out of date.

# Start docker images:

docker run \
  --rm \
  --detach --name cromwell_test_mariadb_database_23306 \
  --env MYSQL_ROOT_PASSWORD=private \
  --env MYSQL_USER=cromwell \
  --env MYSQL_PASSWORD=test \
  --env MYSQL_DATABASE=cromwell_test \
  --publish 23306:3306 \
  --volume ${PWD}/src/ci/docker-compose/mariadb-conf.d:/etc/mysql/conf.d \
  mariadb:5.5

RC=1
while [ "${RC}" != "0" ]
do
  sleep 10
  mysql \
    --protocol=tcp --host=localhost --port=23306 \
    --user=cromwell --password=test \
    --execute='DROP DATABASE IF EXISTS cromwell_test; CREATE DATABASE cromwell_test;'
  RC=$?
done

docker run \
  --rm \
  --detach --name cromwell_test_mariadb_database_33306 \
  --env MYSQL_ROOT_PASSWORD=private \
  --env MYSQL_USER=cromwell \
  --env MYSQL_PASSWORD=test \
  --env MYSQL_DATABASE=cromwell_test \
  --publish 33306:3306 \
  --volume ${PWD}/src/ci/docker-compose/mariadb-conf.d:/etc/mysql/conf.d \
  mariadb:latest
sleep 10

RC=1
while [ "${RC}" != "0" ]
do
  sleep 10
  mysql \
    --protocol=tcp --host=localhost --port=33306 \
    --user=cromwell --password=test \
    --execute='DROP DATABASE IF EXISTS cromwell_test; CREATE DATABASE cromwell_test;'
  RC=$?
done

docker run \
  --rm \
  --detach --name cromwell_test_mysql_database_3306 \
  --env MYSQL_ROOT_PASSWORD=private \
  --env MYSQL_USER=cromwell \
  --env MYSQL_PASSWORD=test \
  --env MYSQL_DATABASE=cromwell_test \
  --publish 3306:3306 \
  --volume ${PWD}/src/ci/docker-compose/mysql-conf.d:/etc/mysql/conf.d \
  mysql:5.6

RC=1
while [ "${RC}" != "0" ]
do
  sleep 10
  mysql \
    --protocol=tcp --host=localhost --port=3306 \
    --user=cromwell --password=test \
    --execute='DROP DATABASE IF EXISTS cromwell_test; CREATE DATABASE cromwell_test;'
  RC=$?
done

docker run \
  --rm \
  --detach --name cromwell_test_mysql_database_13306 \
  --env MYSQL_ROOT_PASSWORD=private \
  --env MYSQL_USER=cromwell \
  --env MYSQL_PASSWORD=test \
  --env MYSQL_DATABASE=cromwell_test \
  --publish 13306:3306 \
  --volume ${PWD}/src/ci/docker-compose/mysql-conf.d:/etc/mysql/conf.d \
  mysql:latest
sleep 10

RC=1
while [ "${RC}" != "0" ]
do
  sleep 10
  mysql \
    --protocol=tcp --host=localhost --port=13306 \
    --user=cromwell --password=test \
    --execute='DROP DATABASE IF EXISTS cromwell_test; CREATE DATABASE cromwell_test;'
  RC=$?
done

docker run \
  --rm \
  --detach --name cromwell_test_postgresql_database_5432 \
  --env POSTGRES_USER=cromwell \
  --env POSTGRES_PASSWORD=test \
  --env POSTGRES_DB=cromwell_test \
  --publish 5432:5432 \
  --volume ${PWD}/src/ci/docker-compose/postgresql-initdb.d:/docker-entrypoint-initdb.d \
  postgres:9.6

RC=1
while [ "${RC}" != "0" ]
do
  sleep 10
  PGPASSWORD=test psql \
    --host=localhost --port=5432 --username=cromwell \
    postgres <<< 'drop database if exists cromwell_test; create database cromwell_test;'
  RC=$?
done

docker run \
  --rm \
  --detach --name cromwell_test_postgresql_database_15432 \
  --env POSTGRES_USER=cromwell \
  --env POSTGRES_PASSWORD=test \
  --env POSTGRES_DB=cromwell_test \
  --publish 15432:5432 \
  --volume ${PWD}/src/ci/docker-compose/postgresql-initdb.d:/docker-entrypoint-initdb.d \
  postgres:latest

RC=1
while [ "${RC}" != "0" ]
do
  sleep 10
  PGPASSWORD=test psql \
    --host=localhost --port=5432 --username=cromwell \
    postgres <<< 'drop database if exists cromwell_test; create database cromwell_test;'
  RC=$?
done

