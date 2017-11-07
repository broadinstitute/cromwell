#!/usr/bin/env bash

set -e

echo "BUILD_TYPE='$BUILD_TYPE'"
echo "TRAVIS_BRANCH='$TRAVIS_BRANCH'"
echo "TRAVIS_PULL_REQUEST='$TRAVIS_PULL_REQUEST'"

publish() {
  docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"
  sbt $@ +publish dockerBuildAndPush
}

if [ "$BUILD_TYPE" == "sbt" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; then

    if [ "$TRAVIS_BRANCH" == "develop" ]; then
        # Publish images for both the "cromwell develop branch" and the "cromwell dev environment".
        CROMWELL_DOCKER_TAGS=develop,dev publish

    elif [[ "$TRAVIS_BRANCH" =~ ^[0-9\.]+_hotfix$ ]]; then
        publish -Dproject.isSnapshot=false
    fi
fi
