import Dependencies.engineDependencies
import Dependencies.coreDependencies
import Merging.customMergeStrategy
import Testing._
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtrelease.ReleasePlugin._

object Settings {
  val engineVersion = "0.20"

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
    test in assembly     := {},
    logLevel in assembly := Level.Info
  )

  val commonSettings = releaseSettings ++ testSettings ++ assemblySettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := "2.11.7",
    resolvers ++= commonResolvers,
    scalacOptions ++= compilerSettings,
    parallelExecution := false
  )

  val coreSettings = List(
    name := "cromwell-core",
    version := engineVersion,
    libraryDependencies ++= coreDependencies,
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar"
  ) ++ commonSettings

  val engineSettings = List(
    name := "cromwell-engine",
    version := engineVersion,
    libraryDependencies ++= engineDependencies,
    assemblyMergeStrategy in assembly := customMergeStrategy,
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar"
  ) ++ commonSettings

  val backendSettings = List(
    name := "cromwell-backend",
    version := "0.1",
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar"
  ) ++ commonSettings

  val rootSettings = List(
    name := "cromwell",
    version := engineVersion,
    assemblyJarName in assembly := name.value + "-" + version.value + ".jar",
    packageOptions in assembly += Package.ManifestAttributes("Premain-Class" -> "org.aspectj.weaver.loadtime.Agent"),
    assemblyMergeStrategy in assembly := customMergeStrategy
  ) ++ commonSettings
}
