#!/usr/bin/env bash

set -e

echo "BUILD_TYPE='$BUILD_TYPE'"
echo "TRAVIS_BRANCH='$TRAVIS_BRANCH'"
echo "TRAVIS_PULL_REQUEST='$TRAVIS_PULL_REQUEST'"

if [ "$BUILD_TYPE" == "sbt" ] && [ "$TRAVIS_BRANCH" == "develop" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; then
    sbt 'set test in Test := {}' publish
fi
