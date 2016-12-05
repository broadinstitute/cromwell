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
    "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/",
    "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/"
  )

  /*
    The reason why -Xmax-classfile-name is set is because this will fail
    to build on Docker otherwise.  The reason why it's 200 is because it
    fails if the value is too close to 256 (even 254 fails).  For more info:

    https://github.com/sbt/sbt-assembly/issues/69
    https://github.com/scala/pickling/issues/10

    Other fancy flags from

    http://blog.threatstack.com/useful-scalac-options-for-better-scala-development-part-1

    and

    https://tpolecat.github.io/2014/04/11/scalac-flags.html

  */
  val compilerSettings = List(
    "-Xlint",
    "-feature",
    "-Xmax-classfile-name", "200",
    "-target:jvm-1.8",
    "-encoding", "UTF-8",
    "-unchecked",
    "-deprecation",
    "-Xfuture",
    "-Yno-adapted-args",
    "-Ywarn-dead-code",
    "-Ywarn-numeric-widen",
    "-Ywarn-value-discard",
    "-Ywarn-unused",
    "-Ywarn-unused-import",
    "-Xfatal-warnings"
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

  val commonSettings = ReleasePlugin.projectSettings ++ testSettings ++ assemblySettings ++
    dockerSettings ++ cromwellVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := "2.11.8",
    resolvers ++= commonResolvers,
    scalacOptions ++= compilerSettings,
    scalacOptions in (Compile, doc) ++= docSettings,
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

  val backendSettings = List(
    name := "cromwell-backend"
  ) ++ commonSettings

  val sfsBackendSettings = List(
    name := "cromwell-sfs-backend"
  ) ++ commonSettings

  val htCondorBackendSettings = List(
    name := "cromwell-htcondor-backend",
    libraryDependencies ++= htCondorBackendDependencies
  ) ++ commonSettings

  val sparkBackendSettings = List(
    name := "cromwell-spark-backend",
    libraryDependencies ++= sparkBackendDependencies
  ) ++ commonSettings

  val jesBackendSettings = List(
    name := "cromwell-jes-backend"
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
