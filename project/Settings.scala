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

  /*
    The reason why -Xmax-classfile-name is set is because this will fail
    to build on Docker otherwise.  The reason why it's 200 is because it
    fails if the value is too close to 256 (even 254 fails).  For more info:

    https://github.com/sbt/sbt-assembly/issues/69
    https://github.com/scala/pickling/issues/10

    Other fancy flags from https://tpolecat.github.io/2017/04/25/scalac-flags.html.

    The following isn't used (yet), and in general is an exercise in pain for 2.12 with Cromwell.
    It'd certainly be nice to have, but params causes a world of hurt. Interested parties are encouraged
    to take a stab at it.

    "-Ywarn-unused:params"              // Warn if a value parameter is unused.
  */
  val compilerSettings = List(
    "-explaintypes",
    "-feature",
    "-Xmax-classfile-name", "200",
    "-target:jvm-1.8",
    "-encoding", "UTF-8",
    "-unchecked",
    "-deprecation",
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

  val docSettings = List(
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    "-no-link-warnings"
  )

  lazy val assemblySettings = Seq(
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar",
    aggregate in assembly := false,
    test in assembly := {},
    logLevel in assembly := Level.Info,
    assemblyMergeStrategy in assembly := customMergeStrategy
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

  val ScalaVersion = "2.12.2"
  val commonSettings = ReleasePlugin.projectSettings ++ testSettings ++ assemblySettings ++
    dockerSettings ++ cromwellVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := ScalaVersion,
    resolvers ++= commonResolvers,
    scalacOptions ++= (compilerSettings ++ consoleHostileSettings),
    scalacOptions in (Compile, doc) ++= docSettings,
    scalacOptions in (Compile, console) := compilerSettings,
    parallelExecution := false
  )

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
    scalacOptions ++= (compilerSettings ++ consoleHostileSettings),
    scalacOptions in (Compile, doc) ++= docSettings,
    scalacOptions in (Compile, console) := compilerSettings,
    resolvers ++= commonResolvers
  ) ++ ReleasePlugin.projectSettings ++ testSettings ++ assemblySettings ++
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

  val rootSettings = List(
    name := "cromwell",
    libraryDependencies ++= rootDependencies
  ) ++ commonSettings

}
