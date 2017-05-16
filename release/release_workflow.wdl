task do_release {
   # Repo to release
   String repo
   # Current version being released
   String releaseV
   # Next version
   String nextV
   # Command that will update the appropriate file for the current release
   String updateVersionCommand
   
   # Commands that will update previously released/published dependencies in this repo
   Array[String] dependencyCommands = []
   
   # When true, nothing will be pushed to github, allows for some level of testing
   Boolean dryRun = false
   
   # Can be swapped out to try this on a fork
   String organization
    
   command {
     set -e
     set -x 
     
     # Clone repo and checkout develop
     git clone https://github.com/${organization}/${repo}.git
     cd ${repo}
     git checkout develop
     git pull --rebase

     # Expect the version number on develop to be the version TO BE RELEASED
     echo "Releasing ${organization}/${repo} ${releaseV}${true=" - This a dry run, push commands won't be executed" false = "" dryRun}"
     
     echo "Updating dependencies"
     ${sep='\n' dependencyCommands}
     
     # Make sure tests pass
     sbt update
     JAVA_OPTS=-XX:MaxMetaspaceSize=512m sbt test
     
     git add .
     # If there is nothing to commit, git commit will return 1 which will fail the script.
     # This ensures we only commit if build.sbt was effectively updated
     git diff-index --quiet HEAD || git commit -m "Update ${repo} version to ${releaseV}"
       
     # wdl4s needs a scala docs update
     if [ ${repo} == "wdl4s" ]; then
     
       # Generate new scaladoc
       sbt 'set scalacOptions in (Compile, doc) := List("-skip-packages", "better")' doc
       git checkout gh-pages
       mv target/scala-2.11/api ${releaseV}
       git add ${releaseV}
       
       # Update latest pointer
       git rm --ignore-unmatch latest
       ln -s ${releaseV} latest
       git add latest
       
       git diff-index --quiet HEAD || git commit -m "Update Scaladoc"
       git push ${true="--dry-run" false="" dryRun} origin gh-pages
       
       # Update badges on README
       git checkout develop
       curl -o scaladoc.png https://img.shields.io/badge/scaladoc-${releaseV}-blue.png
       curl -o version.png https://img.shields.io/badge/version-${releaseV}-blue.png
       
       git add scaladoc.png
       git add version.png
       
       git diff-index --quiet HEAD || git commit -m "Update README badges"
       git push ${true="--dry-run" false="" dryRun} origin develop
     fi
     
     # Merge develop into master
     git checkout master
     git pull --rebase
     git merge develop
     
     # Pin centaur for cromwell
     if [ ${repo} == "cromwell" ]; then
        centaurDevelopHEAD=$(git ls-remote git://github.com/${organization}/centaur.git | grep refs/heads/develop | cut -f 1)
        sed -i '' s/CENTAUR_BRANCH=.*/CENTAUR_BRANCH="$centaurDevelopHEAD"/g .travis.yml
        git add .travis.yml
        git commit -m "Pin release to centaur branch"
     fi 
     
     # Tag the release
     git tag ${releaseV}
     
     # Push master and push the tags
     git push ${true="--dry-run" false="" dryRun} origin master
     git push ${true="--dry-run" false="" dryRun} --tags
     
     # Create and push the hotfix branch
     git checkout -b ${releaseV}_hotfix
     git push origin ${releaseV}_hotfix
     
     # Assemble jar for cromwell
     if [ ${repo} == "cromwell" ]; then
        sbt -Dproject.version=${releaseV} -Dproject.isSnapshot=false assembly
     fi  
     
     # Update develop to point to next release version
     git checkout develop
     ${updateVersionCommand}
     git add .
     git diff-index --quiet HEAD || git commit -m "Update ${repo} version from ${releaseV} to ${nextV}"
     git push ${true="--dry-run" false="" dryRun} origin develop
     
     pwd > executionDir.txt
   }

   output {
     String version = releaseV
     String executionDir = read_string(repo + "/executionDir.txt")
   }
}

task wait_for_artifactory {
    String repo
    String version
    
    command <<<
        checkIfPresent() {
            isPresent=$(curl -s --head https://artifactory.broadinstitute.org/artifactory/simple/libs-release-local/org/broadinstitute/${repo}/${version}/ | head -n 1 | grep -q "HTTP/1.[01] [23]..")
        }
        
        elapsedTime=0
        checkIfPresent
        
        # Allow 20 minutes for the file to appear in artifactory
        while [ $? -ne 0 ] && [ $elapsedTime -lt 1200 ]; do
            sleep 10;
            let "elapsedTime+=10"
            checkIfPresent
        done
        
        exit $?
    >>>
    
    output {
        String publishedVersion = version
    }
}

task create_update_dependency_command {
    String dependencyName
    String newVersion
    String dependencyFilePath = "build.sbt"
    
    command {
        echo "sed -i '' \"s/${dependencyName}[[:space:]]=.*/${dependencyName} = \\\"${newVersion}\\\"/g\" ${dependencyFilePath}"
    }
    
    output {
      String updateCommand = read_string(stdout())
    }
}

task versionPrep {
    String organization
    String repo
    String file
    String regexPrefix
    String updateCommandTemplate
    
    String bash_rematch = "{BASH_REMATCH[1]}"
    command <<<
        curl -o versionFile https://raw.githubusercontent.com/${organization}/${repo}/develop/${file}
        regex="${regexPrefix}\"(([0-9]+\.)?([0-9]+))\""
        
        if [[ $(cat versionFile) =~ $regex ]]
        then
            version="$${bash_rematch}"
            echo $version > version
            echo $version | perl -ne 'if (/^([0-9]+\.)?([0-9]+)$/) { $incr = $2 + 1; print "$1$incr\n" }' > nextVersion
        else
            exit 1
        fi
    >>>
    
    output {
        String version = read_string("version")
        String nextVersion = read_string("nextVersion")
        String updateCommand = sub(updateCommandTemplate, "<<VERSION>>", nextVersion)
    }
}

task makeGithubRelease {
    String githubToken
    String organization
    File cromwellJar
    Int oldVersion
    Int newVersion
    
    command <<<
        set -e
        set -x
        
        # download changelog from master
        curl https://raw.githubusercontent.com/${organization}/cromwell/master/CHANGELOG.md -o CHANGELOG.md
        
        # Extract the latest piece of the changelog corresponding to this release
        # head remove the last line, next sed escapes all ", and last sed/tr replaces all new lines with \n so it can be used as a JSON string
        BODY=$(cat CHANGELOG.md | sed -n '/## ${newVersion}/,/## ${oldVersion}/p' | head -n -1 | sed -e 's/"/\\"/g' | sed 's/$/\\n/' | tr -d '\n')
        
        # Build the json body for the POST release
        API_JSON="{\"tag_name\": \"${newVersion}\",\"name\": \"${newVersion}\",\"body\": \"$BODY\",\"draft\": true,\"prerelease\": false}"
        
        # POST the release as a draft
        curl --data "$API_JSON" https://api.github.com/repos/${organization}/cromwell/releases?access_token=${githubToken} -o release_response
        
        # parse the response to get the release id and the asset upload url
        RELEASE_ID=$(cat release_response | python -c "import sys, json; print json.load(sys.stdin)['id']")
        UPLOAD_URL=$(cat release_response | python -c "import sys, json; print json.load(sys.stdin)['upload_url']")
        # Maybe update this when we have basename !
        UPLOAD_URL=$(echo $UPLOAD_URL | sed 's/{.*}/?name=cromwell-${newVersion}.jar/g')
        
        # Upload the cromwell jar as an asset
        curl -X POST --data-binary @${cromwellJar} -H "Authorization: token ${githubToken}" -H "Content-Type: application/octet-stream" $UPLOAD_URL
        
        # Publish the draft
        curl -X PATCH -d '{"draft": false}' https://api.github.com/repos/horneth/cromwell/releases/$RELEASE_ID?access_token=${githubToken}
    >>>
    runtime {
        docker: "python:2.7"
    }
}

workflow release_cromwell {
  String githubToken
  String organization
    
  Pair[String, String] lenthallAsDependency = ("lenthallV", waitForLenthall.publishedVersion) 
  Pair[String, String] wdl4sAsDependency = ("wdl4sV", waitForWdl4s.publishedVersion)
  
  Array[Pair[String, String]] wdl4sDependencies = [lenthallAsDependency]
  Array[Pair[String, String]] wdltoolDependencies = [wdl4sAsDependency]
  Array[Pair[String, String]] cromwellDependencies = [lenthallAsDependency, wdl4sAsDependency]
  
  # Regex to find the line setting the current version
  String dependencyRegexPrefix = "git\\.baseVersion[[:space:]]:=[[:space:]]"
  # Template command to update the version
  String dependencyTemplate = "sed -i '' \"s/git\\.baseVersion[[:space:]]:=.*/git.baseVersion := \\\"<<VERSION>>\\\",/g\" build.sbt"
  
  String cromwellTemplate = "sed -i '' \"s/cromwellVersion[[:space:]]=.*/cromwellVersion = \\\"<<VERSION>>\\\"/g\" project/Version.scala"
  String cromwellRegexPrefix = "cromwellVersion[[:space:]]=[[:space:]]"
  
  # Prepare releases by finding out the current version, next version, and update version command
  call versionPrep as lenthallPrep { input: 
    organization = organization,
    repo = "lenthall",
    file = "build.sbt",
    regexPrefix = dependencyRegexPrefix,
    updateCommandTemplate = dependencyTemplate
  }
  
  call versionPrep as wdl4sPrep { input: 
    organization = organization,
    repo = "wdl4s",
    file = "build.sbt",
    regexPrefix = dependencyRegexPrefix,
    updateCommandTemplate = dependencyTemplate
  }
  
  call versionPrep as wdltoolPrep { input: 
    organization = organization,
    repo = "wdltool",
    file = "build.sbt",
    regexPrefix = dependencyRegexPrefix,
    updateCommandTemplate = dependencyTemplate
  }
  
  call versionPrep as cromwellPrep { input: 
    organization = organization,
    repo = "cromwell",
    file = "project/Version.scala",
    regexPrefix = cromwellRegexPrefix,
    updateCommandTemplate = cromwellTemplate
  }
  
  # Release calls  
  call do_release as release_lenthall { input: 
        organization = organization, 
        repo = "lenthall", 
        releaseV = lenthallPrep.version,
        nextV = lenthallPrep.nextVersion,
        updateVersionCommand = lenthallPrep.updateCommand,
       }
       
  call do_release as release_wdl4s { input: 
        organization = organization, 
        repo = "wdl4s", 
        releaseV = wdl4sPrep.version,
        nextV = wdl4sPrep.nextVersion,
        updateVersionCommand = wdl4sPrep.updateCommand,
        dependencyCommands = wdl4sDependencyCommands.updateCommand
       }
       
  call do_release as release_wdltool { input:
         organization = organization, 
        repo = "wdltool", 
        releaseV = wdltoolPrep.version,
        nextV = wdltoolPrep.nextVersion,
        updateVersionCommand = wdltoolPrep.updateCommand,
        dependencyCommands = wdltoolDependencyCommands
       }  
       
  call do_release as release_cromwell { input: 
        organization = organization, 
        repo = "cromwell", 
        releaseV = cromwellPrep.version,
        nextV = cromwellPrep.nextVersion,
        updateVersionCommand = cromwellPrep.updateCommand,
        dependencyCommands = cromwellDependencyCommands.updateCommand
       }
  
  call wait_for_artifactory as waitForLenthall { input: repo = "lenthall_2.11", version = release_lenthall.version }
  call wait_for_artifactory as waitForWdl4s { input: repo = "wdl4s_2.11", version = release_wdl4s.version }
  
  # Generates commands to update wdl4s dependencies
  scatter(wdl4sDependency in wdl4sDependencies) {
    String depName = wdl4sDependency.left
    String versionName = wdl4sDependency.right
    
    call create_update_dependency_command as wdl4sDependencyCommands { input: 
       dependencyName = depName,
       newVersion = versionName,
       dependencyFilePath = "build.sbt"
    }
  }
  
  # Generates commands to update wdltool dependencies
  scatter(wdltoolDependency in wdltoolDependencies) {
    String depName = wdltoolDependency.left
    String versionName = wdltoolDependency.right
    
    call create_update_dependency_command as wdltoolDependencyCommands { input: 
       dependencyName = depName,
       newVersion = versionName,
       dependencyFilePath = "build.sbt"
    }
  }
  
  # Generates commands to update cromwell dependencies
  scatter(cromwellDependency in cromwellDependencies) {
    String depName = cromwellDependency.left
    String versionName = cromwellDependency.right
    
    call create_update_dependency_command as cromwellDependencyCommands { input: 
       dependencyName = depName,
       newVersion = versionName,
       dependencyFilePath = "project/Dependencies.scala"
    }
  }
  
  File cromwellJar = release_cromwell.executionDir + "/target/scala-2.11/cromwell-" + cromwellPrep.version + ".jar"
  # Version that was just released
  Int cromwellVersionAsInt = cromwellPrep.version
  # Previous version
  Int cromwellPreviousVersion = cromwellVersionAsInt - 1
  
  call makeGithubRelease { input:
           githubToken = githubToken,
           organization = organization,
           cromwellJar = cromwellJar,
           new_version = cromwellVersionAsInt,
           old_version = cromwellPreviousVersion
  }
  
  output {
    String cromwellJar = release_cromwell.executionDir + "/target/scala-2.11/cromwell-" + cromwellPrep.version + ".jar"
  }
}
