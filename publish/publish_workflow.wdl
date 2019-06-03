# Releases the next version of Cromwell.
# Only handles strictly-increasing versions.
# The current version is retrieved from the GitHub releases page.
# If `majorRelease` is `false`, then `nextVersion = currentVersion + 0.1`.
# If `majorRelease` is `true`, then `nextVersion = floor(currentVersion) + 1`.

version 1.0

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
            ~{if (useEmailForName) then "" else "--user '~{githubName}'"} \

        echo 'Clone repo and checkout develop'
        git clone https://github.com/~{organization}/cromwell.git --branch develop cromwell
        cd cromwell

        echo 'Merge develop into master'
        git checkout --track origin/master
        git merge --no-edit develop

        echo 'Tag the release'
        git tag --message=~{releaseVersion} ~{releaseVersion}

        echo 'Assemble jars for cromwell and womtool'
        sbt -Dproject.version=~{releaseVersion} -Dproject.isSnapshot=false server/assembly womtool/assembly

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
        File cromwellJar = "cromwell/server/target/scala-2.12/cromwell-~{releaseVersion}.jar"
        File womtoolJar = "cromwell/womtool/target/scala-2.12/womtool-~{releaseVersion}.jar"
    }

    runtime {
        docker: publishDocker
        cpu: 4
        memory: "8GB"
    }
}

task doMinorRelease {
    meta {
        description: "simply builds jars from the hotfix branch, tags the branch, and push the tag to github"
    }

    input {
        File githubTokenFile
        String githubUser
        String githubName
        String githubEmail
        String organization
        String releaseVersion
        String publishDocker
    }

    Float releaseVersionAsFloat = releaseVersion
    Int majorReleaseNumber = floor(releaseVersionAsFloat)
    String hotfixBranchName = "~{majorReleaseNumber}_hotfix"
    Boolean useEmailForName = githubName == "null"

    command {
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        echo 'Setup Git'
        /cromwell-publish/git-setup.sh \
            --tokenFile '~{githubTokenFile}' \
            --user '~{githubUser}' \
            --email '~{githubEmail}' \
            ~{if (useEmailForName) then "" else "--user '~{githubName}'"} \

        echo 'Clone repo and checkout hotfix branch'
        git clone https://github.com/~{organization}/cromwell.git --branch ~{hotfixBranchName} cromwell
        cd cromwell

        echo 'Tag the release'
        git tag --message=~{releaseVersion} ~{releaseVersion}

        echo 'Assemble jars for cromwell and womtool'
        sbt -Dproject.version=~{releaseVersion} -Dproject.isSnapshot=false server/assembly womtool/assembly

        echo 'Push the tags'
        git push origin ~{releaseVersion}
    }

    output {
        File cromwellJar = "cromwell/server/target/scala-2.12/cromwell-~{releaseVersion}.jar"
        File womtoolJar = "cromwell/womtool/target/scala-2.12/womtool-~{releaseVersion}.jar"
    }

    runtime {
        docker: publishDocker
        cpu: 4
        memory: "8GB"
    }
}

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
            --fail --silent \
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
            --fail --silent \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --head \
            https://api.github.com/rate_limit \
        | grep -F 'X-OAuth-Scopes:' > scopesHeader.txt

        # All of these scopes are required
        required_scopes=(repo)
        # Plus one scopes of these are required
        email_scopes=(user user:email)

        echo 'Verify that all required scopes are present'
        correct_scopes=true
        for scope in "${required_scopes[@]}"; do
            if ! grep -F " ${scope}" scopesHeader.txt; then
                correct_scopes=false
            fi
        done

        if [[ "${correct_scopes}" == true ]]; then
            # Verify that at least one of the email scopes is present
            correct_scopes=false
            for scope in "${email_scopes[@]}"; do
                if grep -F " ${scope}" scopesHeader.txt; then
                    correct_scopes=true
                fi
            done
        fi

        if [[ "${correct_scopes}" != true ]]; then
             echo "Token must have the scopes: (${required_scopes[*]}), plus one of (${email_scopes[*]})." >&2
             echo "Instead found $(cat scopesHeader.txt || echo "no scopes!")" >&2
             exit 1
        fi

        echo 'Get the GitHub user'
        curl \
            --fail --silent \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            https://api.github.com/user \
        > githubUser.json

        jq --raw-output --exit-status .login < githubUser.json > githubUser.txt

        echo 'Get the list of email addresses registered in GitHub'
        curl \
            --fail --silent \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            https://api.github.com/user/emails \
        > githubEmails.json

        echo 'Get the "Broad" email address, otherwise just use the primary email address'
        jq --raw-output --exit-status 'map(select(.email | contains("broad"))) | .[0] .email' \
        < githubEmails.json > githubEmail.txt \
        || jq --raw-output --exit-status 'map(select(.primary)) | .[0] .email' \
        < githubEmails.json > githubEmail.txt

        echo "Get the user's name, that might be 'null'"
        curl --fail --silent https://api.github.com/users/"$(cat githubUser.txt)" > githubUser.json

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
        String currentReleaseVersion =
            if (majorRelease) then floor(previousReleaseVersionAsFloat) + 1 else previousReleaseVersionAsFloat + 0.1
    }

    runtime {
        docker: publishDocker
    }
}

task draftGithubRelease {
    input {
        File githubTokenFile
        String organization
        Boolean majorRelease
        String previousReleaseVersion
        String currentReleaseVersion
        String publishDocker
    }

    Float previousReleaseVersionAsFloat = previousReleaseVersion
    Int previousMajorReleaseNumber = floor(previousReleaseVersionAsFloat)
    Float currentReleaseVersionAsFloat = currentReleaseVersion
    Int currentMajorReleaseNumber = floor(currentReleaseVersionAsFloat)
    String hotfixBranchName = "~{currentMajorReleaseNumber}_hotfix"
    String changelogPreviousVersion = if (majorRelease) then previousMajorReleaseNumber else previousReleaseVersion
    String changelogBranchName = if (majorRelease) then "develop" else hotfixBranchName

    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        # download changelog from ~{changelogBranchName}
        curl \
            --fail --silent \
            --output CHANGELOG.md \
            https://raw.githubusercontent.com/~{organization}/cromwell/~{changelogBranchName}/CHANGELOG.md

        echo 'Extract the latest piece of the changelog corresponding to this release'

        # sed removes the last line
        sed -n '/## ~{currentReleaseVersion}/,/## ~{changelogPreviousVersion}/p' CHANGELOG.md \
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
            --fail --silent \
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

    String releaseRoot = "https://github.com/~{organization}/cromwell/releases/download/~{releaseVersion}/"

    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        echo 'Upload the cromwell jar as an asset'
        curl \
            --fail --silent \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --header "Content-Type: application/octet-stream" \
            --data-binary @~{cromwellJar} \
            "~{cromwellUploadUrl}"

        echo 'Upload the womtool jar as an asset'
        curl \
            --fail --silent \
            --header "Authorization: token $(cat ~{githubTokenFile})" \
            --header "Content-Type: application/octet-stream" \
            --data-binary @~{womtoolJar} \
            "~{womtoolUploadUrl}"

        echo 'Publish the draft'
        curl \
            --fail --silent \
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

task releaseHomebrew {
    input {
        File githubTokenFile
        String githubUser
        String githubName
        String githubEmail
        String organization
        String releaseVersion
        String cromwellReleaseUrl
        String womtoolReleaseUrl
        File cromwellJar
        File womtoolJar
        String publishDocker
    }

    String branchName = "cromwell-~{releaseVersion}"
    Boolean useEmailForName = githubName == "null"

    # This is to allow for testing. If the organization is not broadinstitute, meaning we're doing a test,
    # instead of creating the PR in the official homebrew repo make a PR in the test repo.
    String headBranch = if (organization == "broadinstitute") then "broadinstitute:~{branchName}" else branchName
    String baseBranch = "master"
    String pullRepo = if (organization == "broadinstitute") then "Homebrew" else organization

    meta {
        doc: "https://docs.brew.sh/How-To-Open-a-Homebrew-Pull-Request"
    }

    # 'brew bump-formula-pr' seems very promising and could simplify a lot of this, however it's unclear if it would be
    # able to update the womtool version too.
    command <<<
        # Do not use `set -x` or it will print the GitHub token!
        set -euo pipefail

        echo 'Setup Git'
        /cromwell-publish/git-setup.sh \
            --tokenFile '~{githubTokenFile}' \
            --user '~{githubUser}' \
            --email '~{githubEmail}' \
            ~{if (useEmailForName) then "" else "--user '~{githubName}'"} \

        echo 'Clone the homebrew fork'
        git clone https://github.com/~{organization}/homebrew-core.git --depth=100
        cd homebrew-core

        echo 'Add the original homebrew repo as a remote'
        # See https://help.github.com/articles/syncing-a-fork/
        git remote add upstream https://github.com/Homebrew/homebrew-core
        echo 'Get the master branch from the original homebrew up to date'
        git fetch upstream
        echo 'Checkout our local master branch'
        git checkout master
        echo 'Merge upstream homebrew to local master'
        git merge --no-edit upstream/master

        echo 'Create branch for release'

        git checkout -b ~{branchName} master

        echo '########################################################'
        echo '##### Update url and hash for cromwell and womtool #####'
        echo '########################################################'

        echo 'Find the line with the cromwell url. Also grab the next one that we assume will be the sha.'
        echo 'Push both in an array'
        cromwell_array=()
        while IFS='' read -r line; do cromwell_array+=("${line}"); done \
        < <( \
            grep \
                -h -A 1 \
                'https://github.com/broadinstitute/cromwell/releases/download/.*/cromwell-.*\.jar' \
                Formula/cromwell.rb \
            | awk -F '"' '{print $2}' \
        )
        # Put the url and sha into separate variables.
        cromwell_previous_url="${cromwell_array[0]}"
        cromwell_previous_sha="${cromwell_array[1]}"

        echo 'Same for womtool'
        womtool_array=()
        while IFS='' read -r line; do womtool_array+=("${line}"); done \
        < <( \
            grep \
                -h -A 1 \
                'https://github.com/broadinstitute/cromwell/releases/download/.*/womtool-.*\.jar' \
                Formula/cromwell.rb \
            | awk -F '"' '{print $2}' \
        )
        womtool_previous_url="${womtool_array[0]}"
        womtool_previous_sha="${womtool_array[1]}"

        echo 'Compute the SHA of the newly released jars'
        cromwell_release_sha=$(shasum -a 256 ~{cromwellJar} | cut -c 1-64)
        womtool_release_sha=$(shasum -a 256 ~{womtoolJar} | cut -c 1-64)

        echo 'Update cromwell url'
        sed -i -e "s;${cromwell_previous_url};~{cromwellReleaseUrl};g" Formula/cromwell.rb
        echo 'Update womtool url'
        sed -i -e "s;${womtool_previous_url};~{womtoolReleaseUrl};g" Formula/cromwell.rb
        echo 'Update cromwell sha'
        sed -i -e "s;${cromwell_previous_sha};${cromwell_release_sha};g" Formula/cromwell.rb
        echo 'Update womtool sha'
        sed -i -e "s;${womtool_previous_sha};${womtool_release_sha};g" Formula/cromwell.rb

        echo '########################################################'
        echo '#####   Verify install works and create PR if so   #####'
        echo '########################################################'

        set +e
        brew uninstall --force cromwell
        brew install --build-from-source Formula/cromwell.rb
        install_exit_status=$?
        brew test Formula/cromwell.rb
        test_exit_status=$?
        brew audit --strict Formula/cromwell.rb
        audit_exit_status=$?
        set -e

        if [[ "${install_exit_status}" -eq 0 && "${test_exit_status}" -eq 0 && "${audit_exit_status}" -eq 0 ]]; then
            echo "Creating Homebrew PR"
            echo 'Download a template for the homebrew PR'
            curl \
                --fail --silent \
                --output template.md \
                https://raw.githubusercontent.com/~{pullRepo}/homebrew-core/master/.github/PULL_REQUEST_TEMPLATE.md

            echo 'Check the boxes in the template. Whitelist the boxes to be checked so that if new boxes are added'
            echo 'they are not automatically checked and should be reviewed manually instead.'
            sed \
                -i \
                -e '/guidelines for contributing/s/\[[[:space:]]\]/[x]/' \
                -e "/checked that there aren't other open/s/\[[[:space:]]\]/[x]/" \
                -e "/built your formula locally/s/\[[[:space:]]\]/[x]/" \
                -e "/your test running fine/s/\[[[:space:]]\]/[x]/" \
                -e "/your build pass/s/\[[[:space:]]\]/[x]/" \
                template.md

            jq \
                --null-input --rawfile template template.md '{
                    title: "Cromwell ~{releaseVersion}",
                    head: "~{headBranch}",
                    base: "~{baseBranch}",
                    body: $template
                }' \
            > template.json

            echo 'Add the changes to our new branch'
            git add Formula/cromwell.rb
            git commit --message "cromwell ~{releaseVersion}"

            echo 'Push the branches'
            git push origin master ~{branchName}

            echo 'Create the pull request'
            curl \
                --fail --silent \
                --header "Content-Type: application/json" \
                --header "Authorization: token $(cat ~{githubTokenFile})" \
                --data @template.json \
                https://api.github.com/repos/~{pullRepo}/homebrew-core/pulls
        else
            echo "Homebrew verification failed." 2>&1
            exit 1
        fi
    >>>

    output {
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
        Boolean majorRelease = true
        Boolean publishHomebrew = true
        String publishDocker = "broadinstitute/cromwell-publish:latest"
    }

    parameter_meta {
        githubToken: "Github token to interact with github API"
        organization:
            "Organization on which the release will be performed. Swap out for a test organization for testing"
        majorRelease: "Set to false to do a .X minor release"
        publishHomebrew: "Set to false for most test cases, especially when testing releases to github only"
    }

    call prepGithub { input:
        githubToken = githubToken,
        organization = organization,
        majorRelease = majorRelease,
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
        majorRelease = majorRelease,
        currentReleaseVersion = currentReleaseVersion,
        previousReleaseVersion = cromwellPreviousVersion,
        publishDocker = publishDocker,
    }

    if (draftGithubRelease.done) {
        if (majorRelease) {
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

        if (!majorRelease) {
            call doMinorRelease { input:
                githubTokenFile = githubTokenFile,
                githubUser = prepGithub.githubUser,
                githubName = prepGithub.githubName,
                githubEmail = prepGithub.githubEmail,
                organization = organization,
                releaseVersion = currentReleaseVersion,
                publishDocker = publishDocker,
            }
        }
    }

    File cromwellJar = select_first([doMajorRelease.cromwellJar, doMinorRelease.cromwellJar])
    File womtoolJar = select_first([doMajorRelease.womtoolJar, doMinorRelease.womtoolJar])

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

    if (publishHomebrew) {
        call releaseHomebrew { input:
            githubTokenFile = githubTokenFile,
            githubUser = prepGithub.githubUser,
            githubName = prepGithub.githubName,
            githubEmail = prepGithub.githubEmail,
            organization = organization,
            releaseVersion = currentReleaseVersion,
            cromwellReleaseUrl = publishGithubRelease.cromwellReleaseUrl,
            womtoolReleaseUrl = publishGithubRelease.womtoolReleaseUrl,
            cromwellJar = cromwellJar,
            womtoolJar = womtoolJar,
            publishDocker = publishDocker,
        }
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
