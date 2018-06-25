# Cromwell Release WDL
# Part 1 of the release process.
# See also: https://docs.google.com/document/d/1khCqpOYpkCE95pE4a6VfBIxd2zHPWuMILgpisNKCF6c

task do_release {
   # Repo to release
   String repo = "cromwell"
   # Current version being released
   String releaseV
   # Next version
   String nextV
   # Command that will update the appropriate file for the current release
   String updateVersionCommand

   # Can be swapped out to try this on a fork
   String organization

   command {
     set -e
     set -x 

     # Clone repo and checkout develop
     git clone git@github.com:${organization}/${repo}.git -b develop ${repo}
     cd ${repo}

     # Expect the version number on develop to be the version TO BE RELEASED

     git add .
     # If there is nothing to commit, git commit will return 1 which will fail the script.
     # This ensures we only commit if build.sbt was effectively updated
     git diff-index --quiet HEAD || git commit -m "Update ${repo} version to ${releaseV}"

     # Merge develop into master
     git checkout -t origin/master
     git merge develop --no-edit

     # Make sure tests pass
     sbt update
     sbt test:compile
     # Test with a backup just in case of timeouts
     # If tests repeatedly timeout locally, ensure all tests pass in Travis and then comment out the next line
     sbt test || sbt testQuick

     # Tag the release
     git tag ${releaseV}

     # Push master and push the tags
     git push origin master
     git push --tags

     # Create and push the hotfix branch
     git checkout -b ${releaseV}_hotfix

     git push origin ${releaseV}_hotfix

     # Assemble jar for cromwell
     sbt -Dproject.version=${releaseV} -Dproject.isSnapshot=false assembly

     # Update develop to point to next release version
     git checkout develop
     ${updateVersionCommand}
     git add .
     git diff-index --quiet HEAD || git commit -m "Update ${repo} version from ${releaseV} to ${nextV}"
     git push origin develop

     pwd > executionDir.txt
   }

   output {
     String version = releaseV
     String executionDir = read_string(repo + "/executionDir.txt")
   }
}

task versionPrep {
    String organization
    String updateCommandTemplate
    String repo = "cromwell"
    String file = "project/Version.scala"

    command <<<
      curl -o versionFile https://raw.githubusercontent.com/${organization}/${repo}/develop/${file}
      perl -ne '/\s+cromwellVersion\s*=\s*"(.*)"/ && print $1' versionFile > currentVersion
      echo $((`cat currentVersion` + 1)) > nextVersion
    >>>

    output {
        String nextVersion = read_string("nextVersion")
        String updateCommand = sub(updateCommandTemplate, "<<VERSION>>", nextVersion)
        String currentVersion = read_string("currentVersion")
    }
}

task makeGithubRelease {
    String githubToken
    String organization
    File cromwellJar
    File womtoolJar
    Int oldVersion
    Int newVersion

    command <<<
        set -e
        set -x

        # download changelog from master
        curl https://raw.githubusercontent.com/${organization}/cromwell/master/CHANGELOG.md -o CHANGELOG.md

        # Extract the latest piece of the changelog corresponding to this release
        # head remove the last line, next sed escapes all ", and last sed/tr replaces all new lines with \n so it can be used as a JSON string
        BODY=$(sed -n '/## ${newVersion}/,/## ${oldVersion}/p' CHANGELOG.md | head -n -1 | sed -e 's/"/\\"/g' | sed 's/$/\\n/' | tr -d '\n')

        # Build the json body for the POST release
        API_JSON="{\"tag_name\": \"${newVersion}\",\"name\": \"${newVersion}\",\"body\": \"$BODY\",\"draft\": true,\"prerelease\": false}"

        # POST the release as a draft
        curl --data "$API_JSON" https://api.github.com/repos/${organization}/cromwell/releases?access_token=${githubToken} -o release_response

        # parse the response to get the release id and the asset upload url
        RELEASE_ID=$(python -c "import sys, json; print json.load(sys.stdin)['id']" < release_response)
        UPLOAD_URL=$(python -c "import sys, json; print json.load(sys.stdin)['upload_url']" < release_response)

        CROMWELL_UPLOAD_URL=$(sed 's/{.*}/?name=cromwell-${newVersion}.jar/' <<< "$UPLOAD_URL")
        WOMTOOL_UPLOAD_URL=$(sed 's/{.*}/?name=womtool-${newVersion}.jar/' <<< "$UPLOAD_URL")

        # Upload the cromwell jar as an asset
        curl -X POST --data-binary @${cromwellJar} -H "Authorization: token ${githubToken}" -H "Content-Type: application/octet-stream" "$CROMWELL_UPLOAD_URL"

        # Upload the womtool jar as an asset
        curl -X POST --data-binary @${womtoolJar} -H "Authorization: token ${githubToken}" -H "Content-Type: application/octet-stream" "$WOMTOOL_UPLOAD_URL"

        # Publish the draft
        curl -X PATCH -d '{"draft": false}' https://api.github.com/repos/${organization}/cromwell/releases/"$RELEASE_ID"?access_token=${githubToken}
    >>>
    runtime {
        docker: "python:2.7"
    }
}

workflow release_cromwell {
  String githubToken
  String organization

  String cromwellTemplate = "sed -i '' \"s/cromwellVersion[[:space:]]=.*/cromwellVersion = \\\"<<VERSION>>\\\"/g\" project/Version.scala"

  call versionPrep { input:
    organization = organization,
    updateCommandTemplate = cromwellTemplate
  }

  # Release calls  
  call do_release { input:
        organization = organization, 
        releaseV = versionPrep.currentVersion,
        nextV = versionPrep.nextVersion,
        updateVersionCommand = versionPrep.updateCommand
       }

  File cromwellJar = "${do_release.executionDir}/server/target/scala-2.12/cromwell-${versionPrep.currentVersion}.jar"
  File womtoolJar = "${do_release.executionDir}/womtool/target/scala-2.12/womtool-${versionPrep.currentVersion}.jar"
  Int cromwellVersionAsInt = versionPrep.currentVersion
  # Previous version
  Int cromwellPreviousVersionAsInt = cromwellVersionAsInt - 1

  call makeGithubRelease { input:
           githubToken = githubToken,
           organization = organization,
           cromwellJar = cromwellJar,
           womtoolJar = womtoolJar,
           newVersion = cromwellVersionAsInt,
           oldVersion = cromwellPreviousVersionAsInt
  }

  output {
    File cromwellReleasedJar = cromwellJar
    File womtoolReleasedJar = womtoolJar
  }
}
