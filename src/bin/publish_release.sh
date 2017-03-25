#!/usr/bin/env bash

set -e

echo "TRAVIS_TAG='$TRAVIS_TAG'"

sbt \
    'set test in Test := {}' \
    "set version := \"${TRAVIS_TAG}\"" \
    'set isSnapshot := false' \
    'set resolvers += Resolver.url("bintray-sbt-plugin-releases", url("http://dl.bintray.com/content/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)' \
    'set publishTo := Option("artifactory-publish" at "https://artifactory.broadinstitute.org/artifactory/libs-release-local;build.timestamp=" + new java.util.Date().getTime)' \
    "set credentials += Credentials(\"Artifactory Realm\", \"artifactory.broadinstitute.org\", \"${ARTIFACTORY_USERNAME}\", \"${ARTIFACTORY_PASSWORD}\")" \
    "+ publish" 2>&1 | sed 's/.*set credentials.*/REDACTED LOG LINE/'
