import Dependencies._
import Merging.customMergeStrategy
import Publishing._
import Testing._
import Version._
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.DockerPlugin.autoImport._
import sbtrelease.ReleasePlugin

object Settings {

  val commonResolvers = List(
    Resolver.jcenterRepo,
    "Broad Artifactory Releases" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
    "Broad Artifactory Snapshots" at "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot/"
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
    "-Ywarn-unused:patvars"
  )

  val consoleHostileSettings = List(
    "-Ywarn-unused:imports", // warns about every unused import on every command.
    "-Xfatal-warnings"       // makes those warnings fatal.
  )

  addCompilerPlugin("org.scalamacros" % "paradise" % "2.1.0" cross CrossVersion.full)

  testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI", "-h", "target/test-reports")

  lazy val assemblySettings = Seq(
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar",
    aggregate in assembly := false,
    test in assembly := {},
    logLevel in assembly := Level.Info,
    assemblyMergeStrategy in assembly := customMergeStrategy
  )

  val ScalaVersion = "2.12.4"
  val commonSettings = ReleasePlugin.projectSettings ++ testSettings ++ assemblySettings ++
    dockerSettings ++ cromwellVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := ScalaVersion,
    resolvers ++= commonResolvers,
    parallelExecution := false,
    scalacOptions ++= (CrossVersion.partialVersion(scalaVersion.value) match {
      case Some((2, 12)) =>
        // The default scalacOptions includes console-hostile options.  These options are overridden specifically below
        // for the `console` target.
        baseSettings ++ warningSettings ++ consoleHostileSettings
      case Some((2, 11)) =>
        // Scala 2.11 takes a simplified set of options
        baseSettings
      case wut => throw new NotImplementedError(s"Found unsupported Scala version $wut. wdl4s does not support versions of Scala other than 2.11 or 2.12.")
    }),
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    scalacOptions in(Compile, doc) := (baseSettings ++ List("-no-link-warnings")),
    // No console-hostile options, otherwise the console is effectively unusable.
    // https://github.com/sbt/sbt/issues/1815
    scalacOptions in(Compile, console) := (baseSettings ++ warningSettings),
    crossScalaVersions := List("2.11.11", "2.12.4")
  )

  val docSettings = List(
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    "-no-link-warnings"
  )

  lazy val dockerSettings = Seq(
    imageNames in docker := Seq(
      ImageName(
        namespace = Option("broadinstitute"),
        repository = name.value,
        tag = Option(cromwellVersion)),
      ImageName(
        namespace = Option("broadinstitute"),
        repository = name.value,
        tag = Option(version.value))
    ),
    dockerfile in docker := {
      // The assembly task generates a fat JAR file
      val artifact: File = assembly.value
      val artifactTargetPath = s"/app/${artifact.name}"

      new Dockerfile {
        from("openjdk:8")
        expose(8000)
        add(artifact, artifactTargetPath)
        runRaw(s"ln -s $artifactTargetPath /app/cromwell.jar")

        // If you use the 'exec' form for an entry point, shell processing is not performed and 
        // environment variable substitution does not occur.  Thus we have to /bin/bash here
        // and pass along any subsequent command line arguments
        // See https://docs.docker.com/engine/reference/builder/#/entrypoint
        entryPoint("/bin/bash", "-c", "java ${JAVA_OPTS} -jar /app/cromwell.jar ${CROMWELL_ARGS} ${*}", "--")
      }
    },
    buildOptions in docker := BuildOptions(
      cache = false,
      removeIntermediateContainers = BuildOptions.Remove.Always
    )
  )

  val lenthallSettings = List(
    name := "cromwell-lenthall",
    libraryDependencies ++= lenthallDependencies
  ) ++ commonSettings

  val womSettings = List(
    name := "cromwell-wom",
    libraryDependencies ++= womDependencies
  ) ++ commonSettings

  val wdlSettings = List(
    name := "cromwell-wdl",
    libraryDependencies ++= wdlDependencies
  ) ++ commonSettings

  val cwlSettings = List(
    name := "cromwell-cwl",
    libraryDependencies ++= cwlDependencies
  ) ++ commonSettings

  val coreSettings = List(
    name := "cromwell-core",
    libraryDependencies ++= coreDependencies
  ) ++ commonSettings

  val servicesSettings = List(
    name := "cromwell-services"
  ) ++ commonSettings

  val gcsFileSystemSettings = List(
    name := "cromwell-gcsfilesystem",
    libraryDependencies ++= gcsFileSystemDependencies
  ) ++ commonSettings

  val databaseSqlSettings = List(
    name := "cromwell-database-sql",
    libraryDependencies ++= databaseSqlDependencies
  ) ++ commonSettings

  val databaseMigrationSettings = List(
    name := "cromwell-database-migration",
    libraryDependencies ++= databaseMigrationDependencies
  ) ++ commonSettings

  val cromwellApiClientSettings = List(
    name := "cromwell-api-client",
    libraryDependencies ++= cromwellApiClientDependencies,
    organization := "org.broadinstitute",
    scalaVersion := ScalaVersion,
    scalacOptions in (Compile, doc) ++= docSettings,
    resolvers ++= commonResolvers
  ) ++ commonSettings ++ ReleasePlugin.projectSettings ++ testSettings ++ assemblySettings ++
    cromwellVersionWithGit ++ publishingSettings

  val dockerHashingSettings = List(
    name := "cromwell-docker-hashing"
  ) ++ commonSettings

  val backendSettings = List(
    name := "cromwell-backend"
  ) ++ commonSettings

  val sfsBackendSettings = List(
    name := "cromwell-sfs-backend"
  ) ++ commonSettings

  val tesBackendSettings = List(
    name := "cromwell-tes-backend",
    libraryDependencies ++= tesBackendDependencies
  ) ++ commonSettings

  val sparkBackendSettings = List(
    name := "cromwell-spark-backend",
    libraryDependencies ++= sparkBackendDependencies
  ) ++ commonSettings

  val jesBackendSettings = List(
    name := "cromwell-jes-backend",
    libraryDependencies ++= jesBackendDependencies
  ) ++ commonSettings

  val engineSettings = List(
    name := "cromwell-engine",
    libraryDependencies ++= engineDependencies
  ) ++ commonSettings ++ versionConfCompileSettings

  val cloudSupportSettings = List(
    name := "cromwell-cloud-support",
    libraryDependencies ++= gcsFileSystemDependencies
  ) ++ commonSettings

  val rootSettings = List(
    name := "cromwell",
    libraryDependencies ++= rootDependencies
  ) ++ commonSettings

}
