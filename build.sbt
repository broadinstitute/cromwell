import sbtassembly.MergeStrategy
import sbtrelease.ReleasePlugin._

name := "cromwell"

version := "0.13"

organization := "org.broadinstitute"

scalaVersion := "2.11.7"

val lenthallV = "0.13-eb2dd2d-SNAPSHOT"

val sprayV = "1.3.2"

val DowngradedSprayV = "1.3.1"

val akkaV = "2.3.12"

val googleClientApiV = "1.20.0"

resolvers ++= Seq(
  "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/",
  "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/")

libraryDependencies ++= Seq(
  "org.broadinstitute" %% "lenthall" % lenthallV,
  "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
  "org.joda" % "joda-convert" % "1.8.1",
  "org.webjars" % "swagger-ui" % "2.1.1",
  "io.spray" %% "spray-can" % sprayV,
  "io.spray" %% "spray-routing" % sprayV,
  "io.spray" %% "spray-client" % sprayV,
  "io.spray" %% "spray-http" % sprayV,
  "io.spray" %% "spray-json" % DowngradedSprayV,
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "com.typesafe.akka" %% "akka-slf4j" % akkaV,
  "commons-codec" % "commons-codec" % "1.10",
  "commons-io" % "commons-io" % "2.4",
  "org.apache.commons" % "commons-lang3" % "3.4",
  "ch.qos.logback" % "logback-classic" % "1.1.3",
  "ch.qos.logback" % "logback-access" % "1.1.3",
  "org.codehaus.janino" % "janino" % "2.7.8",
  "com.typesafe.slick" %% "slick" % "3.0.3",
  "com.zaxxer" % "HikariCP" % "2.3.3",
  "org.hsqldb" % "hsqldb" % "2.3.2",
  "com.google.gcloud" % "gcloud-java" % "latest.integration",
  "com.google.api-client" % "google-api-client-java6" % googleClientApiV,
  "com.google.api-client" % "google-api-client-jackson2" % googleClientApiV,
  "com.google.oauth-client" % "google-oauth-client" % googleClientApiV,
  "mysql" % "mysql-connector-java" % "5.1.36",
  "org.scalaz" % "scalaz-core_2.11" % "7.1.3",
  //---------- Test libraries -------------------//
  "io.spray" %% "spray-testkit" % sprayV % Test,
  "org.scalatest" %% "scalatest" % "2.2.5" % Test,
  "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
  "org.liquibase" % "liquibase-core" % "3.3.5" % Test,
  "org.yaml" % "snakeyaml" % "1.16" % Test
)

releaseSettings

shellPrompt := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value)}

assemblyJarName in assembly := "cromwell-" + version.value + ".jar"

logLevel in assembly := Level.Info

val customMergeStrategy: String => MergeStrategy = {
  case x if Assembly.isConfigFile(x) =>
    MergeStrategy.concat
  case PathList(ps@_*) if Assembly.isReadme(ps.last) || Assembly.isLicenseFile(ps.last) =>
    MergeStrategy.rename
  case PathList("META-INF", path@_*) =>
    path map {
      _.toLowerCase
    } match {
      case ("manifest.mf" :: Nil) | ("index.list" :: Nil) | ("dependencies" :: Nil) =>
        MergeStrategy.discard
      case ps@(x :: xs) if ps.last.endsWith(".sf") || ps.last.endsWith(".dsa") =>
        MergeStrategy.discard
      case "plexus" :: xs =>
        MergeStrategy.discard
      case "spring.tooling" :: xs =>
        MergeStrategy.discard
      case "services" :: xs =>
        MergeStrategy.filterDistinctLines
      case ("spring.schemas" :: Nil) | ("spring.handlers" :: Nil) =>
        MergeStrategy.filterDistinctLines
      case _ => MergeStrategy.deduplicate
    }
  case "asm-license.txt" | "overview.html" | "cobertura.properties" =>
    MergeStrategy.discard
  case _ => MergeStrategy.deduplicate
}

assemblyMergeStrategy in assembly := customMergeStrategy

// The reason why -Xmax-classfile-name is set is because this will fail
// to build on Docker otherwise.  The reason why it's 200 is because it
// fails if the value is too close to 256 (even 254 fails).  For more info:
//
// https://github.com/sbt/sbt-assembly/issues/69
// https://github.com/scala/pickling/issues/10
scalacOptions ++= Seq("-deprecation", "-unchecked", "-feature", "-Xmax-classfile-name", "200")

lazy val DockerTest = config("docker") extend Test

lazy val NoDockerTest = config("nodocker") extend Test

// NOTE: The following block may cause problems with IntelliJ IDEA
// by creating multiple test configurations.
// May need to comment out when importing the project.
lazy val root = (project in file("."))
  .configs(DockerTest).configs(NoDockerTest)
  .settings(inConfig(DockerTest)(Defaults.testTasks): _*)
  .settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)

testOptions in DockerTest += Tests.Argument("-n", "DockerTest")

testOptions in NoDockerTest += Tests.Argument("-l", "DockerTest")

test in assembly := {}

parallelExecution := false

