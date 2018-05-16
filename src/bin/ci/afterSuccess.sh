#!/usr/bin/env bash

set -e

echo "BUILD_TYPE='$BUILD_TYPE'"
echo "TRAVIS_BRANCH='$TRAVIS_BRANCH'"
echo "TRAVIS_EVENT_TYPE='$TRAVIS_EVENT_TYPE'"

publish() {
    # ensure we are not logging while accessing docker
    set +x

    docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"
    sbt $@ +publish dockerBuildAndPush
}

# Some CI environments want to know when new docker images are published. They do not currently poll dockerhub but do
# poll github. To help those environments, signal that a new set of docker images has been published to dockerhub by
# updating a well known branch in github.
push_publish_complete() {
    # ensure we are not logging while accessing vault
    set +x

    # Login to vault to access secrets
    docker run --rm \
        -v $HOME:/root:rw \
        broadinstitute/dsde-toolbox \
        vault auth "$JES_TOKEN" < /dev/null > /dev/null && echo vault auth success

    # Render secrets
    docker run --rm \
        -v $HOME:/root:rw \
        -v $PWD/src/bin/travis/resources:/working \
        -v $PWD:/output \
        -e ENVIRONMENT=not_used \
        -e INPUT_PATH=/working \
        -e OUT_PATH=/output \
        broadinstitute/dsde-toolbox render-templates.sh

    github_private_deploy_key="$(pwd)/github_private_deploy_key"

    # Loosely adapted from https://github.com/broadinstitute/workbench-libs/blob/435a932/scripts/version_update.sh
    mkdir publish_complete
    (
        cd publish_complete

        git init
        git config core.sshCommand "ssh -i $github_private_deploy_key -F /dev/null"
        git config user.email "travis@travis-ci.org"
        git config user.name "Travis CI"

        git_repo="git@github.com:broadinstitute/cromwell.git"
        git_branch="${TRAVIS_BRANCH}_publish_complete"
        git_remote="publish_complete"
        git_message="publish complete [skip ci]"

        git remote add "$git_remote" "$git_repo"
        git checkout -b "$git_branch"
        git commit --allow-empty -m "$git_message"
        git push -f "$git_remote" "$git_branch"
    )
}

if [ "$BUILD_TYPE" == "sbt" ] && [ "$TRAVIS_EVENT_TYPE" == "push" ]; then

    if [ "$TRAVIS_BRANCH" == "develop" ]; then
        # Publish images for both the "cromwell develop branch" and the "cromwell dev environment".
        CROMWELL_DOCKER_TAGS=develop,dev publish
        push_publish_complete

    elif [[ "$TRAVIS_BRANCH" =~ ^[0-9\.]+_hotfix$ ]]; then
        publish -Dproject.isSnapshot=false
    fi
fi
