#!/usr/bin/env bash

set -e
set -x

docker pull ubuntu:latest

ENABLE_COVERAGE=true sbt \
  -Dbackend.providers.Local.config.filesystems.local.localization.0=copy \
  +clean +nointegration:test

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
