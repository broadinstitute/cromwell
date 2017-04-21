#Thoughts--
#Still need to manually check that the jars are being created in broad artifactory
#Need to manually investigate if tests aren't being written
#Need to manualy update the release portion of Cromwell
#Need to make sure that the changelog.md file doesn't get weirdly altered when merging to master for Cromwell


task basicRelease {
   Array[String] dependency
   String repo
   String releaseV
   String baseDir

   command {
     cd ${baseDir}/${repo}
     git checkout develop
     git pull

     currentV=`grep "git.baseVersion " build.sbt | cut -d "=" -f 2 | tr -d \"\,`
     echo "Version of ${repo} on the develop branch is $currentV"

     if [ ${releaseV} -ne $version ]; then
        echo "Need to update the ${repo} version on the develop branch to ${releaseV}" >&2
        sed -i -e 's/'"$currentV"'/'"${releaseV}"'/' build.sbt
        git add build.sbt
        git commit -m "Update ${repo} version from $currentV to ${releaseV}"
        if [ ${repo} -eq "wdl4s" ]; then #updating docs
           sbt 'set scalacOptions in (Compile, doc) := List("-skip-packages", "better")' doc
           git checkout gh-pages
           mv target/scala-2.11/api ${releaseV}
           git add ${releaseV}
           git commit -m "API Docs"
           #git push origin gh-pages

           git checkout develop
           sed -i -e 's/'"$currentV"'/'"${releaseV}"'/' README.md
           git add README.md
           git commit -m "Update ${repo} API docs references"
        fi
        #git push origin develop
     fi

    #git checkout master
    #git merge --squash origin/develop
    #git tag ${releaseV}
    #git push origin master
    #git push --tags
   }

   output{
     String updatedVersion = releaseV
   }
}

task wdltoolRelease {
String baseDir
String repo
String releaseV
String wdl4sV

  command {
    cd ${baseDir}/${repo}
        git checkout develop
        git pull

        currentV=`grep "git.baseVersion " build.sbt | cut -d "=" -f2 | tr -d \"\,`
        echo "Version of ${repo} on the develop branch is $version"

        if [ ${releaseV} -ne $version ];then
           echo "Need to update the ${repo} version on the develop branch to ${releaseV}" >&2
           sed -i -e 's/'"$currentV"'/'"${releaseV}"'/' build.sbt

           #updating the Wdl4s version
           currentWdl4sV=`grep "wdl4s" build.sbt | cut -d "%" -f 4 | tr -d \"\,`
           sed -i -e 's/'"$currentWdl4sV"'/'"${wdl4sV}"'/' build.sbt

           git add build.sbt
           git commit -m "Update ${repo} version from $currentV to ${releaseV}"
           #git push origin develop
        fi

    #git checkout master
    #git merge --squash origin/develop
    #git tag ${releaseV}
    #git push origin master
    #git push --tags
  }
  output {
    String updatedVersion = releaseV
  }

}

task cromwellRelease {
  String repo
  String releaseV
  String baseDir
  String wdl4sV
  String lenthallV

  command {
    cd ${baseDir}/${repo}
    git checkout develop
    git pull

    currentV=`grep "cromwellVersion " project/Version.scala | cut -d "=" -f2 | tr -d \"`
    echo "Version of ${repo} on the develop branch is $version"

    if [ ${releaseV} -ne $version ];then
       echo "Need to update the ${repo} version on the develop branch to ${releaseV}" >&2
       sed -i -e 's/'"$currentV"'/'"${releaseV}"'/' project/Version.scala
       git add project/Version.scala

       #updating the Wdl4s version
       currentWdl4sV=`grep "wdl4sV " project/Dependencies.scala | cut -d "=" -f 2 | tr -d \"\,`
       sed -i -e 's/'"$currentWdl4sV"'/'"${wdl4sV}"'/' project/Dependencies.scala

       #updating the Lenthall version
       currentLenthallV=`grep "lenthallV " project/Dependencies.scala | cut -d "=" -f 2 | tr -d \"\,`
       sed -i -e 's/'"$currentLenthallV"'/'"${lenthallV}"'/' project/Dependencies.scala
       git add project/Dependencies.scala

       #git commit -m "Update ${repo} version from $currentV to ${releaseV} and updated dependencies"
    fi

   #git checkout master
   #git merge --squash origin/develop
   #git tag ${releaseV}
   #git push origin master
   #git push --tags
  }

  output{
    String updatedVersion = releaseV
  }
}

task checkTests {
  String baseDir
  String repo

   command {
      cd ${baseDir}/${repo}
      sbt test > results.txt
      passingTests=`grep "All tests passed." results.txt`

      if [[ -z $passingTests ]]
      then
        echo "Need to fix sbt tests for ${repo}" >&2
        echo "Resuts of Sbt tests are: `cat results.txt`" >&2
        exit 1
      fi
   }
   output {
     String results = read_string("results.txt")
   }
}

task updateDevelopVersion {
  String baseDir
  String repo
  String upcomingV
  String releasedV

  command {
   cd ${baseDir}/${repo}
   git checkout develop

   if [ ${repo} -eq "cromwell ]; then

     sed -i -e 's/'"${releasedV}"'/'"${upcomingV}"'/' project/Version.scala
     git add project/Version.scala

   else

     sed -i -e 's/'"${releasedV}"'/'"${upcomingV}"'/' build.sbt
     git add build.sbt

   fi

   git commit -m "Updating version of upcoming release for ${repo}"
   git push origin develop

  }
  output {}
}


workflow release {

String baseDir = "/Users/rmunshi"

Pair[String,String] lenthall = ("lenthall", "0.20")
Pair[String,String] wdl4s = ("wdl4s", "0.8")
Pair[String,String] wdltool = ("wdltool", "0.8")
Pair[String,String] cromwell = ("cromwell", "24")


Array[Pair[String,Pair[String,String]]] releaseRepos = [("lenthall", ("0.20","0.21")),("wdl4s", ("0.8","0.9")),
                                       ("wdltool", ("0.8","0.9")), ("cromwell", ("24","25"))]

  scatter (repos in releaseRepos) {
    call checkTests { input: baseDir = baseDir, repo=repos.left }
  }

  call basicRelease as lenthallRelease { input:dependdency = checkTests.result
                                               repo = lenthall.left,
                                               releaseV = lenthall.right,
                                               baseDir = baseDir }

  call basicRelease as wdl4sRelease { input:dependency = checkTests.result
                                            repo = wdl4s.left,
                                            releaseV = wdl4s.right,
                                            baseDir = baseDir }

  call wdlToolRelease { input:repo = wdltool.left,
                              releaseV = wdltool.right,
                              baseDir = baseDir
                              wdl4sV = wd4sRelease.updatedVersion }

  call cromwellRelease { input: repo = cromwell.left,
                                releaseV = cromwell.right,
                                baseDir = baseDir,
                                wdl4sV = wdl4sRelease.updatedVersion,
                                lenthallV = lenthallRelease.updatedVersion }


  scatter (repos in upcomingReleases) {
    call updateDevelopVersion { input: repo = repos.left,
                                       baseDir = baseDir,
                                       upcomingV = repos.right.right,
                                       releasedV = repos.right.left }
  }

  output {
    #???
  }
}
