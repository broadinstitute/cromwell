version 1.0

task do_major_release {
   input {
      Int releaseVersion
      String organization
   }
   
   parameter_meta {
     releaseVersion: "Current version being released"
     organization: "Organization on which the release will be performed. Can be used for testing."
   }

   Int nextV = releaseVersion + 1

   command {
     set -e
     set -x 

     # Clone repo and checkout develop
     git clone https://github.com/~{organization}/cromwell.git -b develop cromwell
     cd cromwell

     # Merge develop into master
     git checkout -t origin/master
     git merge develop --no-edit

     # Tag the release
     git tag ~{releaseVersion}

     # Push master and push the tags
     git push origin master
     git push origin ~{releaseVersion}

     # Create and push the hotfix branch
     git checkout -b ~{releaseVersion}_hotfix

     git push origin ~{releaseVersion}_hotfix

     # Assemble jar for cromwell
     sbt -Dproject.version=~{releaseVersion} -Dproject.isSnapshot=false assembly

     # Update develop to point to next release version
     git checkout develop
     sed -i '' -e '/cromwellVersion[[:space:]]=/s/[0-9][0-9]*/~{nextV}/' project/Version.scala
     git add .
     git diff-index --quiet HEAD || git commit -m "Update cromwell version from ~{releaseVersion} to ~{nextV}"
     git push origin develop
   }

   output {
     File cromwellJar = "cromwell/server/target/scala-2.12/cromwell-~{releaseVersion}.jar"
     File womtoolJar = "cromwell/womtool/target/scala-2.12/womtool-~{releaseVersion}.jar"
   }
}

# do_minor_release simply builds jars from the hotfix branch, tags the branch and push the tag to github
task do_minor_release {
   input {
      # Current version being released. e.g: 33.1
      String releaseVersion
   
      # Can be swapped out to try this on a fork
      String organization
   }

   Float releaseVersionAsFloat = releaseVersion
   Int majorReleaseNumber = floor(releaseVersionAsFloat)   
   String hotfixBranchName = "~{majorReleaseNumber}_hotfix"

   command {
     set -e
     set -x 

     # Clone repo and checkout hotfix branch
     git clone https://github.com/~{organization}/cromwell.git -b ~{hotfixBranchName} cromwell
     cd cromwell

     # Make sure tests pass
     sbt update
     sbt test:compile
     # Test with a backup just in case of timeouts
     # If tests repeatedly timeout locally, ensure all tests pass in Travis and then comment out the next line
     sbt test || sbt testQuick

     # Tag the release
     git tag ~{releaseVersion}

     # Push the tags
     git push origin ~{releaseVersion}

     # Assemble jar for cromwell
     sbt -Dproject.version=~{releaseVersion} -Dproject.isSnapshot=false assembly
   }

   output {
     File cromwellJar = "cromwell/server/target/scala-2.12/cromwell-~{releaseVersion}.jar"
     File womtoolJar = "cromwell/womtool/target/scala-2.12/womtool-~{releaseVersion}.jar"
   }
}

task versionPrep {
    input {
        String organization
        Boolean majorRelease
    }

    command <<<
      set -e

      which jq || brew install jq
      curl --fail -v -s https://api.github.com/repos/~{organization}/cromwell/releases/latest | jq --raw-output '.tag_name' > version
    >>>
    
    output {
        String previouslyReleasedVersion = read_string("version")
        Float previouslyReleasedVersionAsFloat = read_float("version")
        String currentReleaseVersion = if (majorRelease) then floor(previouslyReleasedVersionAsFloat) + 1 else previouslyReleasedVersionAsFloat + 0.1
    }
}

task draftGithubRelease {
    input {
        String githubToken
        String organization
        String oldVersion
        String newVersion
    }
    command <<<
        set -e
        set -x

        # download changelog from master
        curl --fail -v https://raw.githubusercontent.com/~{organization}/cromwell/master/CHANGELOG.md -o CHANGELOG.md

        # Extract the latest piece of the changelog corresponding to this release
        # head remove the last line, next sed escapes all ", and last sed/tr replaces all new lines with \n so it can be used as a JSON string
        BODY=$(sed -n '/## ~{newVersion}/,/## ~{oldVersion}/p' CHANGELOG.md | head -n -1 | sed -e 's/"/\\"/g' | sed 's/$/\\n/' | tr -d '\n')

        # Build the json body for the POST release
        API_JSON="{\"tag_name\": \"~{newVersion}\",\"name\": \"~{newVersion}\",\"body\": \"$BODY\",\"draft\": true,\"prerelease\": false}"

        # POST the release as a draft
        curl --fail -v -X POST --data "$API_JSON" https://api.github.com/repos/~{organization}/cromwell/releases?access_token=~{githubToken} -o release_response

        # parse the response to get the release id and the asset upload url
        RELEASE_ID=$(python -c "import sys, json; print json.load(sys.stdin)['id']" < release_response)
        UPLOAD_URL=$(python -c "import sys, json; print json.load(sys.stdin)['upload_url']" < release_response)

        echo $RELEASE_ID > release_id.txt
        echo $UPLOAD_URL > upload_url.txt
    >>>
    runtime {
        docker: "python:2.7"
    }
    output {
      String release_id = read_string("release_id.txt")
      String upload_url = sub(read_string("upload_url.txt"), "\\{.*", "")
      Boolean draft_complete = true
    }

}

task publishGithubRelease {
    input {
        String githubToken
        String organization
        File cromwellJar
        File womtoolJar
        String newVersion

        String release_id
        String upload_url
    }

    String cromwellJarName = basename(cromwellJar)
    String womtoolJarName = basename(womtoolJar)

    String cromwell_upload_url = "~{upload_url}?name=~{cromwellJarName}?label=~{cromwellJarName}"
    String womtool_upload_url = "~{upload_url}?name=~{womtoolJarName}?label=~{womtoolJarName}"

    command <<<
        set -e
        set -x

        cp ~{cromwellJar} crom.jar
        cp ~{womtoolJar} womt.jar

        # Upload the cromwell jar as an asset
        curl --fail -v --data-binary @crom.jar -H "Authorization: token ~{githubToken}" -H "Content-Type: application/octet-stream" "~{cromwell_upload_url}"

        # Upload the womtool jar as an asset
        curl --fail -v --data-binary @womt.jar -H "Authorization: token ~{githubToken}" -H "Content-Type: application/octet-stream" "~{womtool_upload_url}"

        # Publish the draft
        curl --fail -v -X PATCH -d '{"draft": false, "tag": "~{newVersion}"}' https://api.github.com/repos/~{organization}/cromwell/releases/"~{release_id}"?access_token=~{githubToken}
    >>>
    runtime {
        docker: "python:2.7"
    }
    output {
        String cromwellReleaseUrl = "https://github.com/~{organization}/cromwell/releases/download/~{newVersion}/~{cromwellJarName}"
        String womtoolReleaseUrl = "https://github.com/~{organization}/cromwell/releases/download/~{newVersion}/~{womtoolJarName}"
    }
}

task releaseHomebrew {
    input {
        String organization
        String githubToken
        
        String releaseVersion
        
        String cromwellReleaseUrl
        String womtoolReleaseUrl
        
        File cromwellJar
        File womtoolJar
    }
    
    String branchName = "cromwell-~{releaseVersion}"
    
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
        set -e
        set -x 

        # Clone the homebrew fork
        git clone https://github.com/~{organization}/homebrew-core.git --depth=100
        cd homebrew-core
        
        # See https://help.github.com/articles/syncing-a-fork/
        # Add the original homebrew repo as a remote
        git remote add upstream https://github.com/Homebrew/homebrew-core
        # Get the master branch from the original homebrew up to date
        git fetch upstream
        # Checkout our local master branch
        git checkout master
        # Merge upstream homebrew to local master
        git merge upstream/master
        
        # Create branch for release

        git checkout -b ~{branchName} master
        
        ########################################################
        ##### Update url and hash for cromwell and womtool #####
        ########################################################

        # Find the line with the cromwell url. Also grab the next one that we assume will be the sha. Push both in an array 
        IFS=$'\r\n' command eval 'CROMWELL_ARRAY=($(grep -oh -A 1 'https://github.com/broadinstitute/cromwell/releases/download/.*/cromwell-.*\.jar' Formula/cromwell.rb))'
        # Put the url and sha into separate variables. Also trim any leading space
        CROMWELL_PREVIOUS_URL=$(echo "${CROMWELL_ARRAY[0]}" | sed -e 's/^[[:space:]]*//')
        CROMWELL_PREVIOUS_SHA=$(echo "${CROMWELL_ARRAY[1]}" | sed -e 's/^[[:space:]]*//')
        
        # Same for womtool 
        IFS=$'\r\n' command eval 'WOMTOOL_ARRAY=($(grep -oh -A 1 'https://github.com/broadinstitute/cromwell/releases/download/.*/womtool-.*\.jar' Formula/cromwell.rb))'
        WOMTOOL_PREVIOUS_URL=$(echo "${WOMTOOL_ARRAY[0]}" | sed -e 's/^[[:space:]]*//')
        WOMTOOL_PREVIOUS_SHA=$(echo "${WOMTOOL_ARRAY[1]}" | sed -e 's/^[[:space:]]*//')

        # Compute the SHA of the newly released jars
        CROMWELL_RELEASE_SHA=$(shasum -a 256 ~{cromwellJar} | cut -c 1-64)
        WOMTOOL_RELEASE_SHA=$(shasum -a 256 ~{womtoolJar} | cut -c 1-64)

        # Update cromwell url
        sed -i '' -e "s;${CROMWELL_PREVIOUS_URL};~{cromwellReleaseUrl};g" Formula/cromwell.rb
        # Update womtool url
        sed -i '' -e "s;${WOMTOOL_PREVIOUS_URL};~{womtoolReleaseUrl};g" Formula/cromwell.rb
        # Update cromwell sha
        sed -i '' -e "s;${CROMWELL_PREVIOUS_SHA};sha256 \"${CROMWELL_RELEASE_SHA}\";" Formula/cromwell.rb
        # Update womtool sha
        sed -i '' -e "s;${WOMTOOL_PREVIOUS_SHA};sha256 \"${WOMTOOL_RELEASE_SHA}\";" Formula/cromwell.rb
        
        ########################################################
        #####   Verify install works and create PR if so   #####
        ########################################################
        
        brew uninstall --force cromwell
        brew install --build-from-source Formula/cromwell.rb
        INSTALL=$?
        brew audit --strict Formula/cromwell.rb
        AUDIT=$?
        
        if [[ "${INSTALL}" -eq 0 && "${AUDIT}" -eq 0 ]]; then
          git commit -a -m "cromwell ~{releaseVersion}"
          git push origin ~{branchName}
          echo "Creating Homebew PR"
          # Download a template for the homebrew PR
          curl -o template.md https://raw.githubusercontent.com/brewsci/homebrew-bio/master/.github/PULL_REQUEST_TEMPLATE.md
          
          # Check the boxes in the template. White list the boxes to be checked so that if new boxes are added they're not 
          # automatically checked and should be reviewed manually instead
          sed -i '' -e '/guidelines for contributing/s/\[[[:space:]]\]/[x]/' \
          -e "/checked that there aren't other open/s/\[[[:space:]]\]/[x]/" \
          -e "/built your formula locally/s/\[[[:space:]]\]/[x]/" \
          -e "/your build pass/s/\[[[:space:]]\]/[x]/" \
          template.md
          
          curl -H "Content-Type: application/json" \
           -H "Authorization: token ~{githubToken}" \
           -d "{ \"title\": \"Cromwell ~{releaseVersion}\", \"head\": \"~{headBranch}\", \"base\": \"~{baseBranch}\", \"body\": \"$(cat template.md | sed -e 's/"/\\"/g' | sed 's/$/\\n/' | tr -d '\n')\" }" \
           https://api.github.com/repos/~{pullRepo}/homebrew-core/pulls
        fi
    >>>
}

workflow release_cromwell {
  input {
    String githubToken
    String organization = "broadinstitute"
    Boolean majorRelease = true
    Boolean publishHomebrew = true
  }
  
  parameter_meta {
    githubToken: "Github token to interact with github API"
    organization: "Organization on which the release will be performed. Swap out for a test organization for testing"
    majorRelease: "Set to false to do a .X minor release"
  }

  call versionPrep { input:
    organization = organization,
    majorRelease = majorRelease
  }

  call draftGithubRelease { input:
      githubToken = githubToken,
      organization = organization,
      newVersion = cromwellVersion,
      oldVersion = cromwellPreviousVersion
  }

  # This is the version before the one being released
  String cromwellPreviousVersion = versionPrep.previouslyReleasedVersion
  # This is the version being released
  String cromwellVersion = versionPrep.currentReleaseVersion

  if (draftGithubRelease.draft_complete) {
      if (majorRelease) {
        call do_major_release { input:
                organization = organization,
                releaseVersion = cromwellVersion
        }
      }

      if (!majorRelease) {
        call do_minor_release { input:
               organization = organization,
               releaseVersion = cromwellVersion
        }
      }
  }

  File cromwellJar = select_first([do_major_release.cromwellJar, do_minor_release.cromwellJar])
  File womtoolJar = select_first([do_major_release.womtoolJar, do_minor_release.womtoolJar])

  call publishGithubRelease { input:
           githubToken = githubToken,
           organization = organization,
           cromwellJar = cromwellJar,
           womtoolJar = womtoolJar,
           newVersion = cromwellVersion,
           release_id = draftGithubRelease.release_id,
           upload_url = draftGithubRelease.upload_url
  }

  if (publishHomebrew) {
      call releaseHomebrew { input:
               organization = organization,
               githubToken = githubToken,
               releaseVersion = cromwellVersion,
               cromwellReleaseUrl = publishGithubRelease.cromwellReleaseUrl,
               womtoolReleaseUrl = publishGithubRelease.womtoolReleaseUrl,
               cromwellJar = cromwellJar,
               womtoolJar = womtoolJar
      }
  }

  output {
    File cromwellReleasedJar = cromwellJar
    File womtoolReleasedJar = womtoolJar
  }
  
  meta {
    title: "Cromwell Release WDL"
    doc: "https://docs.google.com/document/d/1khCqpOYpkCE95pE4a6VfBIxd2zHPWuMILgpisNKCF6c"
  }
}
