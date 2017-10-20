#!/usr/bin/env bash

set -e

echo "TRAVIS_TAG='$TRAVIS_TAG'"

sbt \
    'set test in Test := {}' \
    "set version := \"${TRAVIS_TAG}\"" \
    'set isSnapshot := false' \
    'set resolvers += Resolver.url("bintray-sbt-plugin-releases", url("http://dl.bintray.com/content/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)' \
    'set publishTo := Option("artifactory-publish" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release-local;build.timestamp=" + new java.util.Date().getTime)' \
    "set credentials += Credentials(\"Artifactory Realm\", \"broadinstitute.jfrog.io\", \"${ARTIFACTORY_USERNAME}\", \"${ARTIFACTORY_PASSWORD}\")" \
    "+ publish" 2>&1 | sed 's/.*set credentials.*/REDACTED LOG LINE/'
