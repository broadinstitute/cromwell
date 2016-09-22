#!/usr/bin/env bash

set -e
set -x

sudo apt-get update -qq
sudo apt-get install -qq mysql-server-5.6 mysql-client-5.6 mysql-client-core-5.6
docker pull ubuntu:latest
mysql -u root -e "SET GLOBAL sql_mode = 'STRICT_ALL_TABLES';"
mysql -u root -e "CREATE DATABASE IF NOT EXISTS cromwell_test;"
mysql -u root -e "CREATE USER 'travis'@'localhost' IDENTIFIED BY '';"
mysql -u root -e "GRANT ALL PRIVILEGES ON cromwell_test . * TO 'travis'@'localhost';"

sbt -Dbackend.providers.Local.config.filesystems.local.localization.0=copy clean coverage nointegration:test coverageReport
sbt coverageAggregate
sbt coveralls
