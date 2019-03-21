import sbt._
import Keys._

val scioVersion = "0.7.3"
val beamVersion = "2.11.0"
val scalaMacrosVersion = "2.1.1"

lazy val commonSettings = Defaults.coreDefaultSettings ++ Seq(
  organization := "example",
  // Semantic versioning http://semver.org/
  version := "0.1.0-SNAPSHOT",
  scalaVersion := "2.12.8",
  scalacOptions ++= Seq("-target:jvm-1.8",
                        "-deprecation",
                        "-feature",
                        "-unchecked"),
  javacOptions ++= Seq("-source", "1.8", "-target", "1.8")
)

lazy val paradiseDependency =
  "org.scalamacros" % "paradise" % scalaMacrosVersion cross CrossVersion.full
lazy val macroSettings = Seq(
  libraryDependencies += "org.scala-lang" % "scala-reflect" % scalaVersion.value,
  addCompilerPlugin(paradiseDependency)
)

lazy val root: Project = project
  .in(file("."))
  .settings(commonSettings)
  .settings(macroSettings)
  .settings(
    name := "scio-job",
    scalacOptions += "-Ypartial-unification",
    description := "scio job",
    publish / skip := true,
    libraryDependencies ++= Seq(
      //"com.spotify" %% "scio-core" % scioVersion,
      "org.typelevel" %% "cats-core" % "1.5.0",
      //"com.spotify" %% "scio-test" % scioVersion % Test,
      "org.apache.beam" % "beam-runners-direct-java" % beamVersion,
      "io.circe" %% "circe-core" % "0.10.0",
      "io.circe" %% "circe-parser" % "0.10.0",
      "io.circe" %% "circe-generic" % "0.10.0",
      "io.circe" %% "circe-java8" % "0.10.0",
      "io.circe" %% "circe-shapes" % "0.10.0",
      "mysql" % "mysql-connector-java" % "6.0.5",
      "com.github.gekomad" %% "itto-csv" % "0.1.0",
      "com.google.code.gson" % "gson" % "2.8.5",

"com.google.cloud.sql" % "mysql-socket-factory-connector-j-6" % "1.0.12",
      // optional dataflow runner
      "org.apache.beam" % "beam-runners-google-cloud-dataflow-java" % beamVersion,
      "org.apache.beam" % "beam-sdks-java-io-jdbc" % beamVersion,
      "org.slf4j" % "slf4j-simple" % "1.7.25"
    )
  )
  .enablePlugins(PackPlugin)

lazy val repl: Project = project
  .in(file(".repl"))
  .settings(commonSettings)
  .settings(macroSettings)
  .settings(
    name := "repl",
    description := "Scio REPL for scio job",
    libraryDependencies ++= Seq(
      "com.spotify" %% "scio-repl" % scioVersion
    ),
    Compile / mainClass := Some("com.spotify.scio.repl.ScioShell"),
    publish / skip := true
  )
  .dependsOn(root)
