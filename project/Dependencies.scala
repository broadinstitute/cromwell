import sbt._

object Dependencies {
  val catsV = "1.0.0-MF"

  val sprayJsonV = "1.3.2"
  val circeVersion = "0.8.0"
  val lenthallV = "0.28-9440cc0-SNAP"

  // Internal collections of dependencies

  private val catsDependencies = List(
    "org.typelevel" %% "cats-core" % catsV,
    "com.github.benhutchison" %% "mouse" % "0.10-MF"
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
    "com.lihaoyi" %% "ammonite-ops" % "1.0.1",
    "org.typelevel" %% "cats-effect" % "0.4",
    "org.pegdown" % "pegdown" % "1.6.0" % Test,
    "org.scalactic" %% "scalactic" % "3.0.1",
    "org.scalatest" %% "scalatest" % "3.0.2" % "test",
    "org.scalacheck" %% "scalacheck" % "1.13.4" % "test"

  ) ++ circeDependencies ++ womDependencies

  private val slf4jBindingDependencies = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-classic" % "1.2.3",
    "ch.qos.logback" % "logback-access" % "1.2.3",
    "org.codehaus.janino" % "janino" % "3.0.7"
  )

  val rootDependencies = slf4jBindingDependencies
}
