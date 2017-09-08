#!/usr/bin/env bash

set -e
set -x

docker pull ubuntu:latest

sbt -Dbackend.providers.Local.config.filesystems.local.localization.0=copy clean coverage nointegration:test coverageReport
sbt coverageAggregate
sbt coveralls
