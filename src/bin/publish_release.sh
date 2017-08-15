#!/usr/bin/env bash

set -e

echo "TRAVIS_TAG='$TRAVIS_TAG'"

sbt \
    'set test in Test := {}' \
    "set version in ThisBuild := \"${TRAVIS_TAG}\"" \
    'set isSnapshot in ThisBuild := false' \
    'set resolvers in ThisBuild += Resolver.url("bintray-sbt-plugin-releases", url("http://dl.bintray.com/content/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)' \
    'set publishTo in ThisBuild := Option("artifactory-publish" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release-local;build.timestamp=" + new java.util.Date().getTime)' \
    "set credentials in ThisBuild += Credentials(\"Artifactory Realm\", \"broadinstitute.jfrog.io\", \"${ARTIFACTORY_USERNAME}\", \"${ARTIFACTORY_PASSWORD}\")" \
    "+ publish" 2>&1 | sed 's/.*set credentials.*/REDACTED LOG LINE/'
