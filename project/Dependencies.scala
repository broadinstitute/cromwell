import sbt._

object Dependencies {
  val catsV = "0.9.0"

  val sprayJsonV = "1.3.2"
  val circeVersion = "0.8.0"
  val lenthallV = "0.25"

  // Internal collections of dependencies

  private val catsDependencies = List(
    "org.typelevel" %% "cats" % catsV,
    "com.github.benhutchison" %% "mouse" % "0.9"
  ) map (_
    /*
    Exclude test framework cats-laws and its transitive dependency scalacheck.
    If sbt detects scalacheck, it tries to run it.
    Explicitly excluding the two problematic artifacts instead of including the three (or four?).
    https://github.com/typelevel/cats/tree/v0.7.2#getting-started
    Re "_2.11" and "_2.12", see also: https://github.com/sbt/sbt/issues/1518
     */
    exclude("org.typelevel", "cats-laws_2.11")
    exclude("org.typelevel", "cats-kernel-laws_2.11")
    exclude("org.typelevel", "cats-laws_2.12")
    exclude("org.typelevel", "cats-kernel-laws_2.12")
    )

  val womDependencies = List(
    "org.broadinstitute" %% "lenthall" % lenthallV,
    "com.typesafe.scala-logging" %% "scala-logging" % "3.5.0",
    "io.spray" %% "spray-json" % sprayJsonV,
    "commons-codec" % "commons-codec" % "1.10",
    "commons-io" % "commons-io" % "2.5",
    "org.apache.commons" % "commons-lang3" % "3.4",
    "com.github.pathikrit" %% "better-files" % "2.17.1",
    "org.scalatest" %% "scalatest" % "3.0.2" % "test"
  ) ++ catsDependencies

  val wdlDependencies = List() ++ womDependencies

  private val circeDependencies = List(
    "generic",
    "generic-extras",
    "shapes",
    "refined",
    "literal"
  ).map(m => "io.circe" %% s"circe-$m" % circeVersion)

  val cwlDependencies = List(
    "io.circe" %% "circe-yaml" % "0.6.1",
    "eu.timepit" %% "refined"            % "0.8.2",
    "com.lihaoyi" %% "ammonite-ops" % "1.0.0-RC7" % "test",
    "org.pegdown" % "pegdown" % "1.6.0" % Test
  ) ++ circeDependencies ++ womDependencies

  private val slf4jBindingDependencies = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-classic" % "1.2.3",
    "ch.qos.logback" % "logback-access" % "1.2.3",
    "org.codehaus.janino" % "janino" % "3.0.7"
  )

  val rootDependencies = slf4jBindingDependencies
}
