import Dependencies._
import Merging.customMergeStrategy
import Publishing._
import Testing._
import Version._
import org.apache.ivy.plugins.conflict.NoConflictManager
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.DockerPlugin
import sbtdocker.DockerPlugin.autoImport._
import sbtrelease.ReleasePlugin
import scoverage.ScoverageKeys._

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
    "-target:jvm-1.8",
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
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar",
    test in assembly := {},
    logLevel in assembly := Level.Error,
    assemblyMergeStrategy in assembly := customMergeStrategy.value//,
    //logLevel in assembly :=
      //sys.env.get("ASSEMBLY_LOG_LEVEL").flatMap(Level.apply).getOrElse((logLevel in assembly).value)
  )

  val Scala2_12Version = "2.12.4"
  val ScalaVersion = Scala2_12Version
  val sharedSettings = ReleasePlugin.projectSettings ++
    cromwellVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := ScalaVersion,
    resolvers ++= commonResolvers,
    parallelExecution := false,
    dependencyOverrides ++= cromwellDependencyOverrides,
    scalacOptions ++= (CrossVersion.partialVersion(scalaVersion.value) match {
      case Some((2, 12)) =>
        // The default scalacOptions includes console-hostile options.  These options are overridden specifically below
        // for the `console` target.
        baseSettings ++ warningSettings ++ consoleHostileSettings
      case _ =>
        throw new NotImplementedError(
          s"Found unsupported Scala version '${scalaVersion.value}'." +
            s" ${name.value} does not support versions of Scala other than 2.12.")
    }),
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    scalacOptions in(Compile, doc) ++= baseSettings ++ List("-no-link-warnings"),
    // No console-hostile options, otherwise the console is effectively unusable.
    // https://github.com/sbt/sbt/issues/1815
    scalacOptions in(Compile, console) --= consoleHostileSettings,
    //
    /*
    Only enable coverage for 2.12.

    NOTE: Like below, gave up coming with an SBT setting. Using an environment variable instead.

    Once 2.11 is gone, instead of
      `ENABLE_COVERAGE=true sbt +test coverageReport`
    one can run
      `sbt coverage test coverageReport`
     */
    coverageEnabled := (CrossVersion.partialVersion(scalaVersion.value) match {
      case Some((2, 12)) => sys.env.get("ENABLE_COVERAGE").exists(_.toBoolean)
      case _ =>
        throw new NotImplementedError(
          s"Found unsupported Scala version '${scalaVersion.value}'." +
            s" ${name.value} does not support versions of Scala other than 2.12.")
    }),
    addCompilerPlugin(paradisePlugin cross CrossVersion.full)
  )

  val crossVersionSettings = List(
    crossScalaVersions := List(Scala2_12Version)
  )

  val dockerTags = settingKey[Seq[String]]("The tags for docker builds.")

  lazy val dockerSettings = Seq(
    /*
    NOTE: Gave up fighting with SBT settings. Using an environment variable instead.

    The below "just works", assuming womtool docker image building is also enabled, setting the right image names and
    versions.

    `sbt 'show docker::imageNames'` returns:
      ArrayBuffer(broadinstitute/womtool:30, broadinstitute/womtool:30-c33be41-SNAP)
      ArrayBuffer(broadinstitute/cromwell:30, broadinstitute/cromwell:30-c33be41-SNAP)

    `CROMWELL_DOCKER_TAGS=dev,develop sbt 'show docker::imageNames'` returns:
      ArrayBuffer(broadinstitute/womtool:dev, broadinstitute/womtool:develop)
      ArrayBuffer(broadinstitute/cromwell:dev, broadinstitute/cromwell:develop)
    */
    dockerTags := sys.env.getOrElse("CROMWELL_DOCKER_TAGS", s"$cromwellVersion,${version.value}").split(","),
    imageNames in docker := dockerTags.value map { tag =>
      ImageName(namespace = Option("broadinstitute"), repository = name.value, tag = Option(tag))
    },
    dockerfile in docker := {
      // The assembly task generates a fat JAR file
      val artifact: File = assembly.value
      val artifactTargetPath = s"/app/${artifact.name}"
      val projectName = name.value

      new Dockerfile {
        from("openjdk:8")
        expose(8000)
        add(artifact, artifactTargetPath)
        runRaw(s"ln -s $artifactTargetPath /app/$projectName.jar")

        // If you use the 'exec' form for an entry point, shell processing is not performed and
        // environment variable substitution does not occur.  Thus we have to /bin/bash here
        // and pass along any subsequent command line arguments
        // See https://docs.docker.com/engine/reference/builder/#/entrypoint
        entryPoint(
          "/bin/bash",
          "-c",
          s"java $${JAVA_OPTS} -jar /app/$projectName.jar $${${projectName.toUpperCase}_ARGS} $${*}",
          "--"
        )
      }
    },
    buildOptions in docker := BuildOptions(
      cache = false,
      removeIntermediateContainers = BuildOptions.Remove.Always
    )
  )

  val swaggerUiSettings = List(resourceGenerators in Compile += writeSwaggerUiVersionConf)
  val backendSettings = List(addCompilerPlugin(kindProjectorPlugin))
  val engineSettings = swaggerUiSettings
  val cromiamSettings = swaggerUiSettings

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
                            integrationTests: Boolean = false,
                            crossCompile: Boolean = false): Project = {

      val builders: Seq[Project => Project] = List(
        addTestSettings,
        if (integrationTests) addIntegrationTestSettings else identity,
        _
          .disablePlugins(AssemblyPlugin)
          .settings(if (crossCompile) crossVersionSettings else List.empty)
          .settings(resourceGenerators in Compile += writeProjectVersionConf)
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
                               buildDocker: Boolean = true): Project = {

      val builders: Seq[Project => Project] = List(
        addTestSettings,
        if (buildDocker) _.enablePlugins(DockerPlugin).settings(dockerSettings) else identity,
        _
          .settings(assemblySettings)
          .settings(resourceGenerators in Compile += writeProjectVersionConf)
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
          .settings(GenerateRestApiDocs.generateRestApiDocsSettings)
      )

      buildProject(project, "root", Nil, builders)
    }
  }

}
