#!/usr/bin/env bash

set -e
set -x

docker pull ubuntu:latest

sbt -Dbackend.providers.Local.config.filesystems.local.localization.0=copy clean coverage nointegration:test coverageReport
sbt coverageAggregate
# Disabling broken coveralls on 29_hotfix. develop already using codecov instead.
