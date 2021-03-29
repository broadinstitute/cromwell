import ContinuousIntegration._
import Dependencies._
import GenerateRestApiDocs._
import Merging.customMergeStrategy
import Publishing._
import Testing._
import Version._
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.{DockerPlugin, Instruction, Instructions}

object Settings {

  val commonResolvers = List(
    Resolver.jcenterRepo,
    "Broad Artifactory Releases" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
    "Broad Artifactory Snapshots" at "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot/",
    Resolver.sonatypeRepo("releases")
  )

  /* The reason why -Xmax-classfile-name is set is because this will fail
     to build on Docker otherwise.  The reason why it's 200 is because it
     fails if the value is too close to 256 (even 254 fails).  For more info:

     https://github.com/sbt/sbt-assembly/issues/69
     https://github.com/scala/pickling/issues/10

     Other fancy flags from https://tpolecat.github.io/2017/04/25/scalac-flags.html.

     Per JG's work in Cromwell, the following can't be turned on without causing piles of errors in wdl4s.  Many of the
     constructs that are flagged look suspicious and probably warrant further scrutiny, but no time for that now.

     "-Ywarn-unused:params"              // Warn if a value parameter is unused.
  */
  val baseSettings = List(
    "-unchecked",
    "-deprecation",
    "-feature",
    "-explaintypes",
    "-Xmax-classfile-name", "200",

    // the backend runs bytecode serialization, classfile writing and method-local
    // optimizations (-opt:l:method) in parallel on N threads
    "-Ybackend-parallelism", "3",
    "-Ycache-plugin-class-loader:last-modified",
    "-Ycache-macro-class-loader:last-modified",
    "-encoding", "UTF-8"
  )

  val warningSettings = List(
    "-Xfuture",
    "-Xlint:adapted-args",
    "-Xlint:by-name-right-associative",
    "-Xlint:constant",
    "-Xlint:delayedinit-select",
    "-Xlint:doc-detached",
    "-Xlint:inaccessible",
    "-Xlint:infer-any",
    "-Xlint:missing-interpolator",
    "-Xlint:nullary-override",
    "-Xlint:nullary-unit",
    "-Xlint:option-implicit",
    "-Xlint:package-object-classes",
    "-Xlint:poly-implicit-overload",
    "-Xlint:private-shadow",
    "-Xlint:stars-align",
    "-Xlint:type-parameter-shadow",
    "-Xlint:unsound-match",
    "-Yno-adapted-args",
    "-Ywarn-dead-code",
    "-Ywarn-numeric-widen",
    "-Ywarn-value-discard",
    "-Ywarn-inaccessible",
    "-Ywarn-unused:implicits",
    "-Ywarn-unused:privates",
    "-Ywarn-unused:locals",
    "-Ypartial-unification",
    "-Ywarn-unused:patvars"
  )

  val consoleHostileSettings = List(
    "-Ywarn-unused:imports", // warns about every unused import on every command.
    "-Xfatal-warnings"       // makes those warnings fatal.
  )

  lazy val assemblySettings = Seq(
    assembly / assemblyJarName := name.value + "-" + version.value + ".jar",
    assembly / test := {},
    assembly / assemblyMergeStrategy := customMergeStrategy.value,
    assembly / logLevel :=
      sys.env.get("CROMWELL_SBT_ASSEMBLY_LOG_LEVEL").flatMap(Level.apply).getOrElse((assembly / logLevel).value)
  )

  // 2.12.13 blocked on the release of sbt-scoverage 1.6.2 https://github.com/scoverage/sbt-scoverage/issues/319
  val Scala2_12Version = "2.12.12"
  val ScalaVersion = Scala2_12Version
  val sharedSettings =
    cromwellVersionWithGit ++ artifactorySettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := ScalaVersion,
    resolvers ++= commonResolvers,
    // Don't run tasks in parallel, especially helps in low CPU environments like Travis
    Global / parallelExecution := false,
    Global / concurrentRestrictions ++= List(
      // Don't run any other tasks while running tests, especially helps in low CPU environments like Travis
      Tags.exclusive(Tags.Test),
      // Only run tests on one sub-project at a time, especially helps in low CPU environments like Travis
      Tags.limit(Tags.Test, 1)
    ),
    dependencyOverrides ++= cromwellDependencyOverrides,
    scalacOptions ++= baseSettings ++ warningSettings ++ consoleHostileSettings,
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    Compile / doc / scalacOptions ++= baseSettings ++ List("-no-link-warnings"),
    // No console-hostile options, otherwise the console is effectively unusable.
    // https://github.com/sbt/sbt/issues/1815
    Compile / console / scalacOptions --= consoleHostileSettings,
    addCompilerPlugin(paradisePlugin),
    excludeDependencies ++= List(
      "org.typelevel" % "simulacrum-scalafix-annotations_2.12",
      "org.typelevel" % "simulacrum-scalafix-annotations_2.13"
    )
  )

  /*
      Docker instructions to install Google Cloud SDK image in docker image. It also installs `crcmod` which
      is needed while downloading large files using `gsutil`
      References:
        - https://stackoverflow.com/questions/28372328/how-to-install-the-google-cloud-sdk-in-a-docker-image
        - https://cromwell.readthedocs.io/en/develop/backends/Google/#issues-with-composite-files
        - https://cloud.google.com/storage/docs/gsutil/addlhelp/CRC32CandInstallingcrcmod
   */
  val installGcloudSettings: List[Def.Setting[Seq[Instruction]]] = List(
    dockerCustomSettings := List(
      Instructions.Env("PATH", "$PATH:/usr/local/gcloud/google-cloud-sdk/bin"),
      // instructions to install `crcmod`
      Instructions.Run("apt-get -y update"),
      Instructions.Run("apt-get -y install python3.8"),
      Instructions.Run("apt -y install python3-pip"),
      Instructions.Run("apt-get -y install gcc python-dev python-setuptools"),
      Instructions.Run("pip3 uninstall crcmod"),
      Instructions.Run("pip3 install --no-cache-dir -U crcmod"),
      Instructions.Run("update-alternatives --install /usr/bin/python python /usr/bin/python3 1"),
      Instructions.Env("CLOUDSDK_PYTHON", "python3"),
      // instructions to install Google Cloud SDK
      Instructions.Run("curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz"),
      Instructions.Run("""mkdir -p /usr/local/gcloud \
                         | && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
                         | && /usr/local/gcloud/google-cloud-sdk/install.sh""".stripMargin),
    )
  )

  val swaggerUiSettings = List(Compile / resourceGenerators += writeSwaggerUiVersionConf)
  val backendSettings = List(addCompilerPlugin(kindProjectorPlugin))
  val engineSettings = swaggerUiSettings
  val cromiamSettings = swaggerUiSettings
  val drsLocalizerSettings = installGcloudSettings

  private def buildProject(project: Project,
                           projectName: String,
                           dependencies: Seq[ModuleID],
                           builders: Seq[Project => Project]): Project = {
    val projectSettings = List(
      name := projectName,
      libraryDependencies ++= dependencies
    ) ++ sharedSettings

    builders.foldRight(project.settings(projectSettings))(_ (_))
  }

  // Adds settings to build a library
  implicit class ProjectLibrarySettings(val project: Project) extends AnyVal {
    def withLibrarySettings(libraryName: String,
                            dependencies: Seq[ModuleID] = List.empty,
                            customSettings: Seq[Setting[_]] = List.empty,
                            integrationTests: Boolean = false): Project = {

      val builders: Seq[Project => Project] = List(
        addTestSettings,
        if (integrationTests) addIntegrationTestSettings else identity,
        _
          .disablePlugins(AssemblyPlugin)
          .settings(Compile / resourceGenerators += writeProjectVersionConf)
          .settings(customSettings)
      )

      buildProject(project, libraryName, dependencies, builders)
    }
  }

  // Adds settings to build an executable
  implicit class ProjectExecutableSettings(val project: Project) extends AnyVal {
    def withExecutableSettings(executableName: String,
                               dependencies: Seq[ModuleID] = List.empty,
                               customSettings: Seq[Setting[_]] = List.empty,
                               buildDocker: Boolean = true,
                               pushDocker: Boolean = true): Project = {

      val builders: Seq[Project => Project] = List(
        addTestSettings,
        if (buildDocker) {
          _
            .enablePlugins(DockerPlugin)
            .settings(dockerSettings)
            .settings(dockerPushSettings(pushDocker))
        } else {
          identity
        },
        _
          .settings(assemblySettings)
          .settings(Compile / resourceGenerators += writeProjectVersionConf)
          .settings(customSettings)
      )

      buildProject(project, executableName, dependencies, builders)
    }
  }

  // Adds settings to build the root project
  implicit class ProjectRootSettings(val project: Project) extends AnyVal {
    def withRootSettings(): Project = {

      val builders: Seq[Project => Project] = List(
        addTestSettings,
        _
          .disablePlugins(AssemblyPlugin)
          .settings(publish := {})
          .settings(generateRestApiDocsSettings)
          .settings(ciSettings)
          .settings(rootArtifactorySettings)
      )

      buildProject(project, "root", Nil, builders)
    }

    /**
      * After aggregations have been added to the root project, we can do additional tasks like checking if every
      * sub-project in build.sbt will also be tested by the root-aggregated `sbt test` command.
      */
    def withAggregateSettings(): Project = {
      project.settings(aggregateSettings(project))
    }
  }

}
