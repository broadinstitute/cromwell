#!/usr/bin/env bash

set -e

echo "BUILD_TYPE='$BUILD_TYPE'"
echo "TRAVIS_TAG='$TRAVIS_TAG'"

if [ "$BUILD_TYPE" == "sbt" ]; then
    sbt -Dproject.version="${TRAVIS_TAG}" -Dproject.isSnapshot=false 'set test in Test := {}' publish
fi
