#!/usr/bin/env bash

set -e

echo "TRAVIS_BRANCH='$TRAVIS_BRANCH'"
echo "TRAVIS_PULL_REQUEST='$TRAVIS_PULL_REQUEST'"

# For now, obfuscate SNAPSHOTs from sbt's developers: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then
    if [ "$TRAVIS_BRANCH" == "develop" ]; then
        sbt \
            'set test in Test := {}' \
            'set isSnapshot := false' \
            'set resolvers += Resolver.url("bintray-sbt-plugin-releases", url("http://dl.bintray.com/content/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)' \
            'set publishTo := Option("artifactory-publish" at "https://artifactory.broadinstitute.org/artifactory/libs-release-local;build.timestamp=" + new java.util.Date().getTime)' \
            "set credentials += Credentials(\"Artifactory Realm\", \"artifactory.broadinstitute.org\", \"${ARTIFACTORY_USERNAME}\", \"${ARTIFACTORY_PASSWORD}\")" \
            "+ publish" 2>&1 | sed 's/.*set credentials.*/REDACTED LOG LINE/'

    elif [[ "$TRAVIS_BRANCH" =~ ^[0-9\.]+_hotfix$ ]]; then
        sbt \
            'set test in Test := {}' \
            'set isSnapshot := false' \
            'set git.gitUncommittedChanges := false' \
            'set resolvers += Resolver.url("bintray-sbt-plugin-releases", url("http://dl.bintray.com/content/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)' \
            'set publishTo := Option("artifactory-publish" at "https://artifactory.broadinstitute.org/artifactory/libs-release-local;build.timestamp=" + new java.util.Date().getTime)' \
            "set credentials += Credentials(\"Artifactory Realm\", \"artifactory.broadinstitute.org\", \"${ARTIFACTORY_USERNAME}\", \"${ARTIFACTORY_PASSWORD}\")" \
            "+ publish" 2>&1 | sed 's/.*set credentials.*/REDACTED LOG LINE/'

    fi
fi
