task do_release {
   # Repo to release
   String repo
   # Version number that should be released
   String releaseV
   # Next version number after this released is done (in preparation for next release)
   String nextReleaseV
   # Command that will update the appropriate file for the current release
   String updateVersionForCurrentRelaseCommand
   # Command that will update the appropriate file for the next release
   String updateVersionForNextRelaseCommand
   # A command to optionnally create a fat jar
   String? assembleJarCommand
   
   # Commands that will update previously released/pubslished dependencies in this repo
   Array[String] dependencyCommands = []
   
   # When true, nothing will be pushed to github, allows for some level of testing
   Boolean dryRun = false
   
   # Can be swapped out to try this on a fork
   String organization = "broadinstitute"
    
   command {
     set -e
     set -x 
     
     # Clone repo and checkout develop
     git clone https://github.com/${organization}/${repo}.git
     cd ${repo}
     git checkout develop
     git pull --rebase

     # Extract current version  
     currentV=$(grep "git.baseVersion " build.sbt | cut -d "=" -f 2 | tr -d \"\,[:space:])

     echo "Updating ${repo} to ${releaseV}"
     ${updateVersionForCurrentRelaseCommand}
     
     echo "Update dependencies"
     ${sep='\n' dependencyCommands}
     
     # Make sure it compiles
     sbt update
     sbt compile
     
     git add .
     # If there is nothing to commit, git commit will return 1 which will fail the script.
     # This ensures we only commit if build.sbt was effectively updated
     git diff-index --quiet HEAD || git commit -m "Update ${repo} version from $currentV to ${releaseV}"
       
     # wdl4s needs a scala docs update
     if [ ${repo} == "wdl4s" ]; then
       sbt 'set scalacOptions in (Compile, doc) := List("-skip-packages", "better")' doc
       git checkout gh-pages
       mv target/scala-2.11/api ${releaseV}
       git add ${releaseV}
       git commit -m "API Docs"
       git push ${true="--dry-run" false="" dryRun} origin gh-pages
       git checkout develop
       sed -i '' s/$currentV/${releaseV}/g README.md
       git add README.md
       git diff-index --quiet HEAD || git commit -m "Update ${repo} API docs references"
     fi
     
     # Merge develop into master
     git checkout master
     git pull --rebase
     git merge develop
     
     # Tag the release
     git tag ${releaseV}
     
     # Push master and push the tags
     git push ${true="--dry-run" false="" dryRun} origin master
     git push ${true="--dry-run" false="" dryRun} --tags
     
     # Assemble jar if required
     ${default = "" assembleJarCommand}
     pwd > executionDir.txt
     
     # Update build.sbt on develop to point to next release version
     git checkout develop
     ${updateVersionForNextRelaseCommand}
     git add .
     git diff-index --quiet HEAD || git commit -m "Update ${repo} version from ${releaseV} to ${nextReleaseV}"
     git push ${true="--dry-run" false="" dryRun} origin develop
     
     # Create and push the hotfix branch
     git checkout -b ${releaseV}_hotfix ${releaseV}
     git push origin ${releaseV}_hotfix
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
            isPresent=$(curl -s --head https://artifactory.broadinstitute.org/artifactory/simple/libs-release-local/org/broadinstitute/${repo}/${version}/ | head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null)
        }
        
        elapsedTime=0
        checkIfPresent
        
        # Allow 10 minutes for the file to appear in artifactory
        while [ $? -ne 0 ] && [ $elapsedTime -lt 600 ]; do
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

task create_update_version_command {
    String newVersion
    String versionFilePath = "build.sbt"
    
    command {
        echo "sed -i '' \"s/git.baseVersion[[:space:]]:=.*/git.baseVersion := \\\"${newVersion}\\\",/g\" ${versionFilePath}"
    }
    
    output {
      String updateCommand = read_string(stdout())
    }
}

task create_update_version_command_for_cromwell {
    String newVersion
    
    command {
        echo "sed -i '' \"s/cromwellVersion[[:space:]]=.*/cromwellVersion = \\\"${newVersion}\\\"/g\" project/Version.scala"
    }
    
    output {
      String updateCommand = read_string(stdout())
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

workflow release_cromwell {
  # Lenthall
  String lenthallRelease
  String lenthallNext
  Pair[String, String] lenthallAsDependency = ("lenthallV", waitForLenthall.publishedVersion) 
  
  # Wdl4s
  String wdl4sRelease
  String wdl4sNext
  Pair[String, String] wdl4sAsDependency = ("wdl4sV", waitForWdl4s.publishedVersion)
  # Wdl4s depends on lenthall
  Array[Pair[String, String]] wdl4sDependencies = [lenthallAsDependency]
  
  #Wdltool
  String wdltoolRelease
  String wdltoolNext
  # Wdltool depends on wdl4s
  Array[Pair[String, String]] wdltoolDependencies = [wdl4sAsDependency]
  
  # Cromwell
  String cromwellRelease
  String cromwellNext
  # Cromwell depends on lenthall and wdl4s
  Array[Pair[String, String]] cromwellDependencies = [lenthallAsDependency, wdl4sAsDependency]
    
  # Release calls  
  call do_release as release_lenthall { input: 
        repo = "lenthall", 
        releaseV = lenthallRelease,
        nextReleaseV = lenthallNext,
        updateVersionForCurrentRelaseCommand = updateLenthallCurrent.updateCommand,
        updateVersionForNextRelaseCommand = updateLenthallNext.updateCommand
       }
       
  call do_release as release_wdl4s { input: 
        repo = "wdl4s", 
        releaseV = wdl4sRelease,
        nextReleaseV = wdl4sNext,
        updateVersionForCurrentRelaseCommand = updateWdl4sCurrent.updateCommand,
        updateVersionForNextRelaseCommand = updateWdl4sNext.updateCommand,
        dependencyCommands = wdl4sDependencyCommands.updateCommand
       }
       
  call do_release as release_wdltool { input: 
        repo = "wdltool", 
        releaseV = wdl4sRelease,
        nextReleaseV = wdl4sNext,
        updateVersionForCurrentRelaseCommand = updateWdltoolCurrent.updateCommand,
        updateVersionForNextRelaseCommand = updateWdltoolNext.updateCommand,
        dependencyCommands = wdltoolDependencyCommands
       }  
       
  call do_release as release_cromwell { input: 
        repo = "cromwell", 
        releaseV = cromwellRelease,
        nextReleaseV = cromwellNext,
        updateVersionForCurrentRelaseCommand = updateCromwellCurrent.updateCommand,
        updateVersionForNextRelaseCommand = updateCromwellNext.updateCommand,
        dependencyCommands = cromwellDependencyCommands.updateCommand,
        assembleJarCommand = "sbt -Dproject.version=" + cromwellRelease + " -Dproject.isSnapshot=false assembly"
       }
  
  call create_update_version_command as updateLenthallCurrent { input: newVersion = lenthallRelease }
  call create_update_version_command as updateLenthallNext { input: newVersion = lenthallNext }
  call create_update_version_command as updateWdl4sCurrent { input: newVersion = wdl4sRelease }
  call create_update_version_command as updateWdl4sNext { input: newVersion = wdl4sNext }
  call create_update_version_command as updateWdltoolCurrent { input: newVersion = wdltoolRelease }
  call create_update_version_command as updateWdltoolNext { input: newVersion = wdltoolNext }
  
  call create_update_version_command_for_cromwell as updateCromwellCurrent { input: newVersion = cromwellRelease }
  call create_update_version_command_for_cromwell as updateCromwellNext { input: newVersion = cromwellNext }
  
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
  
  output {
    String cromwellJar = release_cromwell.executionDir + "/target/scala-2.11/cromwell-" + cromwellRelease + ".jar"
  }
}
