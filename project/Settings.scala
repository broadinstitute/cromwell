import Dependencies._
import Merging.customMergeStrategy
import Publishing._
import Testing._
import Version._
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtrelease.ReleasePlugin._

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
    //"-deprecation", // TODO: PBE: Re-enable deprecation warnings
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

  val commonSettings = releaseSettings ++ testSettings ++ assemblySettings ++
    cromwellVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := "2.11.7",
    resolvers ++= commonResolvers,
    scalacOptions ++= compilerSettings,
    parallelExecution := false
  )

  val coreSettings = List(
    name := "cromwell-core",
    libraryDependencies ++= coreDependencies
  ) ++ commonSettings

  val servicesSettings = List(
    name := "cromwell-services",
    libraryDependencies ++= serviceDependencies
  ) ++ commonSettings

  val gcsFileSystemSettings = List(
    name := "cromwell-gcsfilesystem",
    libraryDependencies ++= gcsFileSystemDependencies
  ) ++ commonSettings

  val databaseSettings = List(
    name := "cromwell-database",
    libraryDependencies ++= databaseDependencies
  ) ++ commonSettings

  val backendSettings = List(
    name := "cromwell-backend"
  ) ++ commonSettings

  val localBackendSettings = List(
    name := "cromwell-local-backend"
  ) ++ commonSettings

  val htCondorBackendSettings = List(
    name := "cromwell-htcondor-backend",
    libraryDependencies ++= htCondorBackendDependencies
  ) ++ commonSettings

  val sgeBackendSettings = List(
    name := "cromwell-sge-backend"
  ) ++ commonSettings

  val jesBackendSettings = List(
    name := "cromwell-jes-backend"
  ) ++ commonSettings

  val engineSettings = List(
    name := "cromwell-engine",
    libraryDependencies ++= engineDependencies
  ) ++ commonSettings

  val rootSettings = List(
    name := "cromwell"
  ) ++ commonSettings
}
