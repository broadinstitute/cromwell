#!/usr/bin/env bash

set -e
set -o xtrace

# sbt publish publishes libs to Artifactory for the scala version sbt is running as.
cd codegen_java
if [[ "$TRAVIS_PULL_REQUEST" == "false" && "$TRAVIS_BRANCH" == "develop" ]]; then
	sbt --warn -Dproject.isSnapshot=false "+ publishLocal"
else
	sbt --warn -Dproject.isSnapshot=true "+ publishLocal"
fi
