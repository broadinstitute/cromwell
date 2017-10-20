#!/usr/bin/env bash

set -e
set -x

docker pull ubuntu:latest

sbt -Dbackend.providers.Local.config.filesystems.local.localization.0=copy clean coverage nointegration:test coverageReport
sbt "project lenthall" "++ 2.11.11" test
sbt "project wom" "++ 2.11.11" test
sbt "project wdl" "++ 2.11.11" test
sbt "project cwl" "++ 2.11.11" test
sbt coverageAggregate
bash <(curl -s https://codecov.io/bash)
