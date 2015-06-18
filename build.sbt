import sbtassembly.Plugin.AssemblyKeys._
import sbtassembly.Plugin._
import sbtrelease.ReleasePlugin._

name := "cromwell"

organization := "org.broadinstitute"

scalaVersion := "2.11.6"

val sprayV = "1.3.2"
val DowngradedSprayV = "1.3.1"
val akkaV = "2.3.6"

libraryDependencies ++= Seq(
  "com.gettyimages" %% "spray-swagger" % "0.5.1",
  "org.webjars" % "swagger-ui" % "2.0.24",
  "io.spray" %% "spray-can" % sprayV,
  "io.spray" %% "spray-routing" % sprayV,
  "io.spray" %% "spray-client" % sprayV,
  "io.spray" %% "spray-http" % sprayV,
  "io.spray" %% "spray-json" % DowngradedSprayV,
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "commons-codec" % "commons-codec" % "1.10",
  "ch.qos.logback" % "logback-classic" % "1.1.3",
  //---------- Test libraries -------------------//
  "io.spray" %% "spray-testkit" % sprayV % Test,
  "org.scalatest" %% "scalatest" % "2.2.5" % Test,
  "com.typesafe.akka" %% "akka-testkit" % akkaV % Test
)

releaseSettings

shellPrompt := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value)}

jarName in assembly := "cromwell-" + version.value + ".jar"

logLevel in assembly := Level.Info

val customMergeStrategy: String => MergeStrategy = {
  case x if Assembly.isConfigFile(x) =>
    MergeStrategy.concat
  case PathList(ps@_*) if (Assembly.isReadme(ps.last) || Assembly.isLicenseFile(ps.last)) =>
    MergeStrategy.rename
  case PathList("META-INF", xs@_*) =>
    xs map {
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
  case "asm-license.txt" | "overview.html" =>
    MergeStrategy.discard
  case _ => MergeStrategy.deduplicate
}

mergeStrategy in assembly := customMergeStrategy

scalacOptions ++= Seq("-deprecation", "-unchecked", "-feature")

lazy val DockerTest = config("docker") extend(Test)

lazy val NoDockerTest = config("nodocker") extend(Test)

lazy val root = (project in file("."))
  .configs(DockerTest).configs(NoDockerTest)
  .settings(inConfig(DockerTest)(Defaults.testTasks): _*)
  .settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)  

testOptions in DockerTest := Seq(Tests.Argument("-n", "DockerTest"))

testOptions in NoDockerTest := Seq(Tests.Argument("-l", "DockerTest"))
