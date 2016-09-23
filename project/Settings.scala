import Dependencies._
import Merging.customMergeStrategy
import Publishing._
import Testing._
import Version._
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtrelease.ReleasePlugin
import sbtdocker.DockerPlugin.autoImport._

object Settings {

  val commonResolvers = List(
    "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/",
    "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/"
  )

  /*
    The reason why -Xmax-classfile-name is set is because this will fail
    to build on Docker otherwise.  The reason why it's 200 is because it
    fails if the value is too close to 256 (even 254 fails).  For more info:

    https://github.com/sbt/sbt-assembly/issues/69
    https://github.com/scala/pickling/issues/10
  */
  val compilerSettings = List(
    "-deprecation",
    "-unchecked",
    "-feature",
    "-Xmax-classfile-name",
    "200"
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
          tag = Some(s"${version.value}")
        )
    ),
    dockerfile in docker := {
      // The assembly task generates a fat JAR file
      val artifact: File = assembly.value
      val artifactTargetPath = s"/app/${artifact.name}"

      new Dockerfile {
        from("openjdk:8")
        expose(8000)
        add(artifact, artifactTargetPath)
        entryPoint("/bin/bash", "-c", "java -jar " + artifactTargetPath + " $CROMWELL_ARGS")
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
  ) ++ commonSettings

  val rootSettings = List(
    name := "cromwell",
    libraryDependencies ++= rootDependencies
  ) ++ commonSettings

}
