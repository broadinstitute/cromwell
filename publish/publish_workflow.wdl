# Releases the next version of Cromwell.
# Only handles strictly-increasing versions.
# The current version is retrieved from the GitHub releases page.
# `nextVersion = floor(currentVersion) + 1`.

version 1.0

# Outputs: JAR files
# Github: git tag, hotfix branch, update version number, merge `develop` into `master`
task doMajorRelease {
    input {
        File githubTokenFile
        String githubUser
        String githubName
        String githubEmail
        String organization
        Int releaseVersion
        String publishDocker
    }

    parameter_meta {
        releaseVersion: "Current version being released"
        organization: "Organization on which the release will be performed. Can be used for testing."
    }

    Int nextVersion = releaseVersion + 1
    Boolean useEmailForName = githubName == "null"

    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        echo 'Setup Git'
        /cromwell-publish/git-setup.sh \
            --tokenFile '~{githubTokenFile}' \
            --user '~{githubUser}' \
            --email '~{githubEmail}' \
            ~{if (useEmailForName) then "" else "--name '~{githubName}'"} \

        echo 'Clone repo and checkout develop'
        git clone https://github.com/~{organization}/cromwell.git --branch develop cromwell
        cd cromwell

        echo 'Merge develop into master'
        git checkout --track origin/master
        git merge --no-edit develop

        echo 'Tag the release'
        git tag --message=~{releaseVersion} ~{releaseVersion}

        # Use sbt.server.forcestart to workaround https://github.com/sbt/sbt/issues/6101
        echo 'Assemble jars for cromwell and womtool'
        sbt \
        -Dsbt.server.forcestart=true \
        -Dproject.version=~{releaseVersion} \
        -Dproject.isSnapshot=false \
        server/assembly \
        womtool/assembly

        echo 'Smoke test Cromwell and Womtool artifacts'
        cat > hello.wdl <<FIN
        task hello {
          String name

          command {
            echo 'hello \${name}!'
          }
          output {
            File response = stdout()
          }
        }

        workflow test {
          call hello
        }
        FIN

        cat > hello.inputs <<FIN
        {
          "test.hello.name": "world"
        }
        FIN

        java -jar server/target/scala-2.13/cromwell-~{releaseVersion}.jar run --inputs hello.inputs hello.wdl
        java -jar womtool/target/scala-2.13/womtool-~{releaseVersion}.jar validate --inputs hello.inputs hello.wdl

        echo 'Create the hotfix branch'
        git checkout -b ~{releaseVersion}_hotfix

        echo 'Update develop to point to next release version'
        git checkout develop
        sed -i -e '/cromwellVersion[[:space:]]=/s/[0-9][0-9]*/~{nextVersion}/' project/Version.scala
        git add project/Version.scala
        git diff-index --quiet HEAD \
        || git commit --message "Update cromwell version from ~{releaseVersion} to ~{nextVersion}"

        echo 'Push branches and tags'
        git push origin master develop ~{releaseVersion} ~{releaseVersion}_hotfix
    >>>

    output {
        File cromwellJar = "cromwell/server/target/scala-2.13/cromwell-~{releaseVersion}.jar"
        File womtoolJar = "cromwell/womtool/target/scala-2.13/womtool-~{releaseVersion}.jar"
    }

    runtime {
        docker: publishDocker
        cpu: 4
        memory: "8GB"
    }
}

# Outputs: parameters for Github
# Github: get version number, check token scopes, user info
# - Consider simplifying version check, latest by time is probably fine now
task prepGithub {
    input {
        String githubToken
        String organization
        Boolean majorRelease
        String publishDocker
    }

    # Use a File so that the github token is never written to the stdout/stderr
    File githubTokenFile = write_lines([githubToken])

    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        curl \
            --location --fail --silent --show-error \
            https://api.github.com/repos/~{organization}/cromwell/releases \
        > releases.json

        # Latest is the latest release by time, not version number.
        # Instead, grab the first page of releases, then find the max release version from the returned results
        # https://rosettacode.org/wiki/Determine_if_a_string_is_numeric
        jq \
            --raw-output --exit-status '
            def is_numeric:
                true and try tonumber catch false;
                [ .[] | select(.tag_name|is_numeric) ] | max_by(.tag_name|tonumber) | .tag_name
            ' \
        < releases.json > version.txt

        echo 'Verify that the token has the scopes that will be required later'

        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --head \
            https://api.github.com/rate_limit \
        | grep -i -F 'X-OAuth-Scopes:' > scopesHeader.txt || true

        # All of these scopes are required with an exact match
        required_scopes=(read:org repo user:email workflow)

        echo 'Verify that all required scopes are present'
        correct_scopes=true
        for scope in "${required_scopes[@]}"; do
            if ! grep -F " ${scope}" scopesHeader.txt; then
                correct_scopes=false
            fi
        done

        if [[ "${correct_scopes}" != true ]]; then
             echo "Token must have the exact scopes: (${required_scopes[*]})" >&2
             echo "Instead found $(cat scopesHeader.txt || echo "no scopes!")" >&2
             exit 1
        fi

        echo 'Get the GitHub user'
        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            https://api.github.com/user \
        > githubUser.json

        jq --raw-output --exit-status .login < githubUser.json > githubUser.txt

        echo 'Get the list of email addresses registered in GitHub'
        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            https://api.github.com/user/emails \
        > githubEmails.json

        echo 'Get the "Broad" email address, otherwise just use the primary email address'
        jq --raw-output --exit-status 'map(select(.email | contains("broad"))) | .[0] .email' \
        < githubEmails.json > githubEmail.txt \
        || jq --raw-output --exit-status 'map(select(.primary)) | .[0] .email' \
        < githubEmails.json > githubEmail.txt

        echo "Get the user's name, that might be 'null'"
        curl --location --fail --silent --show-error https://api.github.com/users/"$(cat githubUser.txt)" > githubUser.json

        jq --raw-output --exit-status .name < githubUser.json > githubName.txt \
        || echo "null" > githubName.txt
    >>>

    output {
        String previousReleaseVersion = read_string("version.txt")
        Float previousReleaseVersionAsFloat = read_float("version.txt")
        File githubTokenFileOutput = githubTokenFile
        String githubUser = read_string("githubUser.txt")
        String githubEmail = read_string("githubEmail.txt")
        String githubName = read_string("githubName.txt")
        String currentReleaseVersion = floor(previousReleaseVersionAsFloat) + 1
    }

    runtime {
        docker: publishDocker
    }
}

# Outputs: url to upload files, release ID to finalize later
# Github: prepare release notes, create draft release
task draftGithubRelease {
    input {
        File githubTokenFile
        String organization
        String previousReleaseVersion
        String currentReleaseVersion
        String publishDocker
    }

    Float previousReleaseVersionAsFloat = previousReleaseVersion
    Int previousMajorReleaseNumber = floor(previousReleaseVersionAsFloat)
    Float currentReleaseVersionAsFloat = currentReleaseVersion
    Int currentMajorReleaseNumber = floor(currentReleaseVersionAsFloat)
    String hotfixBranchName = "~{currentMajorReleaseNumber}_hotfix"
    String changelogPreviousVersion = previousMajorReleaseNumber
    String changelogBranchName = "develop"

    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        # download changelog from ~{changelogBranchName}
        curl \
            --location --fail --silent --show-error \
            --output CHANGELOG.md \
            https://raw.githubusercontent.com/~{organization}/cromwell/~{changelogBranchName}/CHANGELOG.md

        echo 'Extract the latest piece of the changelog corresponding to this release'

        # sed removes the last line
        sed -n '/## ~{currentReleaseVersion}/,/## [0-9]/p' CHANGELOG.md \
        | sed '$d' \
        > changelog.txt

        echo 'Build the json body for the POST release'
        jq \
            --null-input --rawfile changelog changelog.txt '{
                tag_name: "~{currentReleaseVersion}",
                name: "~{currentReleaseVersion}",
                body: $changelog,
                draft: true,
                prerelease: false
            }' \
        > changelog.json

        echo 'POST the release as a draft'
        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --request POST \
            --data @changelog.json \
            --output releaseResponse.json \
            https://api.github.com/repos/~{organization}/cromwell/releases

        echo 'parse the response to get the release id and the asset upload url'
        jq --raw-output --exit-status .id < releaseResponse.json > releaseId.txt
        jq --raw-output --exit-status .upload_url < releaseResponse.json > uploadUrl.txt
    >>>

    output {
        String releaseId = read_string("releaseId.txt")
        String uploadUrl = sub(read_string("uploadUrl.txt"), "\\{.*", "")
        Boolean done = true
    }

    runtime {
        docker: publishDocker
    }
}

# Outputs: JAR URLs
# Github: upload JARs, publish draft
task publishGithubRelease {
    input {
        File githubTokenFile
        String organization
        String releaseVersion
        String releaseId
        String uploadUrl
        File cromwellJar
        File womtoolJar
        String publishDocker
    }

    String cromwellJarName = basename(cromwellJar)
    String womtoolJarName = basename(womtoolJar)

    String cromwellUploadUrl = "~{uploadUrl}?name=~{cromwellJarName}&label=~{cromwellJarName}"
    String womtoolUploadUrl = "~{uploadUrl}?name=~{womtoolJarName}&label=~{womtoolJarName}"

    String releaseRoot = "https://github.com/~{organization}/cromwell/releases/download/~{releaseVersion}"

    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        echo 'Upload the cromwell jar as an asset'
        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --header "Content-Type: application/octet-stream" \
            --data-binary @~{cromwellJar} \
            "~{cromwellUploadUrl}"

        echo 'Upload the womtool jar as an asset'
        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --header "Content-Type: application/octet-stream" \
            --data-binary @~{womtoolJar} \
            "~{womtoolUploadUrl}"

        echo 'Publish the draft'
        curl \
            --location --fail --silent --show-error \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --request PATCH \
            --data '{"draft": false, "tag": "~{releaseVersion}"}' \
            https://api.github.com/repos/~{organization}/cromwell/releases/"~{releaseId}"
    >>>

    output {
        String cromwellReleaseUrl = "~{releaseRoot}/~{cromwellJarName}"
        String womtoolReleaseUrl = "~{releaseRoot}/~{womtoolJarName}"
        Boolean done = true
    }

    runtime {
        docker: publishDocker
    }
}

workflow publish_workflow {
    input {
        String githubToken
        String organization
        String publishDocker = "broadinstitute/cromwell-publish:latest"
    }

    parameter_meta {
        githubToken: "Github token to interact with github API"
        organization:
            "Organization on which the release will be performed. Swap out for a test organization for testing"
    }

    call prepGithub { input:
        githubToken = githubToken,
        organization = organization,
        publishDocker = publishDocker,
    }

    # Use a File so that the github token is never written to the stdout/stderr
    File githubTokenFile = prepGithub.githubTokenFileOutput
    # This is the version before the one being released
    String cromwellPreviousVersion = prepGithub.previousReleaseVersion
    # This is the version being released
    String currentReleaseVersion = prepGithub.currentReleaseVersion

    call draftGithubRelease { input:
        githubTokenFile = githubTokenFile,
        organization = organization,
        currentReleaseVersion = currentReleaseVersion,
        previousReleaseVersion = cromwellPreviousVersion,
        publishDocker = publishDocker,
    }

    if (draftGithubRelease.done) {
        call doMajorRelease { input:
            githubTokenFile = githubTokenFile,
            githubUser = prepGithub.githubUser,
            githubName = prepGithub.githubName,
            githubEmail = prepGithub.githubEmail,
            organization = organization,
            releaseVersion = currentReleaseVersion,
            publishDocker = publishDocker,
        }
    }

    File cromwellJar = select_first([doMajorRelease.cromwellJar])
    File womtoolJar = select_first([doMajorRelease.womtoolJar])

    call publishGithubRelease { input:
        githubTokenFile = githubTokenFile,
        organization = organization,
        releaseVersion = currentReleaseVersion,
        releaseId = draftGithubRelease.releaseId,
        uploadUrl = draftGithubRelease.uploadUrl,
        cromwellJar = cromwellJar,
        womtoolJar = womtoolJar,
        publishDocker = publishDocker,
    }

    output {
        File cromwellReleasedJar = cromwellJar
        File womtoolReleasedJar = womtoolJar
    }

    meta {
        title: "Cromwell Publish WDL"
        doc: "https://docs.google.com/document/d/1khCqpOYpkCE95pE4a6VfBIxd2zHPWuMILgpisNKCF6c"
    }
}
