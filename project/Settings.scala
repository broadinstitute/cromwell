import ContinuousIntegration._
import Dependencies._
import GenerateRestApiDocs._
import Merging._
import Publishing._
import Testing._
import Version._
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.{DockerPlugin, Instruction, Instructions}

object Settings {

  /*
     Fancy flags from https://tpolecat.github.io/2017/04/25/scalac-flags.html.

     Per JG's work in Cromwell, the following can't be turned on without causing piles of errors in wdl4s.  Many of the
     constructs that are flagged look suspicious and probably warrant further scrutiny, but no time for that now.

     "-Ywarn-unused:params"              // Warn if a value parameter is unused.
  */
  val baseSettings = List(
    "-unchecked",
    "-deprecation",
    "-feature",
    "-explaintypes",
    // the backend runs bytecode serialization, classfile writing and method-local
    // optimizations (-opt:l:method) in parallel on N threads
    "-Ybackend-parallelism", "3",
    "-Ycache-plugin-class-loader:last-modified",
    "-Ycache-macro-class-loader:last-modified",
    "-encoding", "UTF-8",
    "-Ymacro-annotations"
  )

  val warningSettings = List(
    "-Xlint:adapted-args",
    "-Xlint:constant",
    "-Xlint:delayedinit-select",
    "-Xlint:doc-detached",
    "-Xlint:inaccessible",
    "-Xlint:infer-any",
    "-Xlint:missing-interpolator",
    "-Xlint:nullary-unit",
    "-Xlint:option-implicit",
    "-Xlint:package-object-classes",
    "-Xlint:poly-implicit-overload",
    "-Xlint:private-shadow",
    "-Xlint:stars-align",
    "-Xlint:type-parameter-shadow",
    "-Ywarn-dead-code",
    "-Ywarn-numeric-widen",
    "-Ywarn-value-discard",
    "-Ywarn-unused:implicits",
    "-Ywarn-unused:privates",
    "-Ywarn-unused:locals",
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
  )

  val Scala2_13Version = "2.13.9"
  private val ScalaVersion: String = Scala2_13Version
  private val sharedSettings: Seq[Setting[_]] =
    cromwellVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := ScalaVersion,
    resolvers ++= additionalResolvers,
    Global / parallelExecution := true,
    // Seems to cause race conditions in tests, that are not pertinent to what's being tested
    Test / parallelExecution := false,
    Global / concurrentRestrictions ++= List(
      // Don't run any other tasks while running tests
      Tags.exclusive(Tags.Test),
      // Only run tests on one sub-project at a time
      Tags.limit(Tags.Test, 1)
    ),
    dependencyOverrides ++= cromwellDependencyOverrides,
    excludeDependencies ++= cromwellExcludeDependencies,
    scalacOptions ++= baseSettings ++ warningSettings ++ consoleHostileSettings,
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    Compile / doc / scalacOptions ++= baseSettings ++ List("-no-link-warnings"),
    // No console-hostile options, otherwise the console is effectively unusable.
    // https://github.com/sbt/sbt/issues/1815
    Compile / console / scalacOptions --= consoleHostileSettings,
    excludeDependencies ++= List(
      "org.typelevel" % "simulacrum-scalafix-annotations_2.13"
    )
  )

  val pact4sSettings = sharedSettings ++ List(
    libraryDependencies ++= pact4sDependencies,
    /**
      * Invoking pact tests from root project (sbt "project pact" test)
      * will launch tests in a separate JVM context that ensures contracts
      * are written to the pact/target/pacts folder. Otherwise, contracts
      * will be written to the root folder.
      */
    Test / fork := true
  ) ++ assemblySettings

  /*
      Docker instructions to install Google Cloud SDK image in docker image. It also installs `crcmod` which
      is needed while downloading large files using `gsutil`.
      References:
        - https://stackoverflow.com/questions/28372328/how-to-install-the-google-cloud-sdk-in-a-docker-image
        - https://cromwell.readthedocs.io/en/develop/backends/Google/#issues-with-composite-files
        - https://cloud.google.com/storage/docs/gsutil/addlhelp/CRC32CandInstallingcrcmod

      Install `getm` for performant signed URL downloading with integrity checking.
      References:
        - https://github.com/xbrianh/getm
   */
  val installLocalizerSettings: List[Setting[Seq[Instruction]]] = List(
    dockerCustomSettings := List(
      Instructions.Env("PATH", "$PATH:/usr/local/gcloud/google-cloud-sdk/bin"),
      // instructions to install `crcmod`
      Instructions.Run("apt-get -y update"),
      Instructions.Run("apt-get -y install python3.11"),
      Instructions.Run("apt-get -y install python3-pip"),
      Instructions.Run("apt-get -y install wget gcc python3-dev python3-setuptools"),
      Instructions.Run("pip3 uninstall crcmod"),
      Instructions.Run("pip3 install --no-cache-dir -U crcmod"),
      Instructions.Run("update-alternatives --install /usr/bin/python python /usr/bin/python3 1"),
      Instructions.Env("CLOUDSDK_PYTHON", "python3"),
      // instructions to install Google Cloud SDK
      Instructions.Run("wget https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz -O /tmp/google-cloud-sdk.tar.gz"),
      Instructions.Run("""mkdir -p /usr/local/gcloud \
                         | && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
                         | && /usr/local/gcloud/google-cloud-sdk/install.sh""".stripMargin),
      // instructions to install `getm`. Pin to version 0.0.5 as the behaviors of future versions with respect to
      // messages or exit codes may change.
      Instructions.Run("pip3 install getm==0.0.5")
    )
  )

  val swaggerUiSettings = List(Compile / resourceGenerators += writeSwaggerUiVersionConf)
  val backendSettings = List(addCompilerPlugin(kindProjectorPlugin))
  val engineSettings: List[Setting[_]] = swaggerUiSettings
  val cromiamSettings: List[Setting[_]] = swaggerUiSettings
  val drsLocalizerSettings: List[Setting[_]] = installLocalizerSettings

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
          .settings(rootPublishingSettings)
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
