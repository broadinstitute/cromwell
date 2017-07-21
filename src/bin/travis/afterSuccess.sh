#!/usr/bin/env bash

set -e

echo "BUILD_TYPE='$BUILD_TYPE'"
echo "TRAVIS_BRANCH='$TRAVIS_BRANCH'"
echo "TRAVIS_PULL_REQUEST='$TRAVIS_PULL_REQUEST'"

if [ "$BUILD_TYPE" == "sbt" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; then

    if [ "$TRAVIS_BRANCH" == "develop" ]; then
        # Publish images for both the "cromwell develop branch" and the "cromwell dev environment".
        docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"
        sbt \
          'set test in Test := {}' \
          'set imageNames in docker := Seq("develop", "dev").map(tag => ImageName(s"broadinstitute/cromwell:$tag"))' \
          publish \
          dockerBuildAndPush

    elif [[ "$TRAVIS_BRANCH" =~ ^[0-9\.]+_hotfix$ ]]; then
        docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"
        sbt 'set test in Test := {}' -Dproject.isSnapshot=false dockerBuildAndPush

    fi

fi
