name := "CromIam"
organization := "org.broadinstitute"
version := "1.0"

scalaVersion := "2.12.1"

scalacOptions := List(
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

scalacOptions in (Compile, doc) ++= List(
  // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
  "-no-link-warnings"
)

libraryDependencies ++= {
  val akkaV       = "2.4.17"
  val akkaHttpV   = "10.0.6"
  val scalaTestV  = "3.0.1"
  val catsV = "0.9.0"
  val lenthallV = "0.24-c8af383-SNAP"

  val catsDependencies = List(
    "org.typelevel" %% "cats" % catsV
  ) map (_
    /*
    Exclude test framework cats-laws and its transitive dependency scalacheck.
    If sbt detects scalacheck, it tries to run it.
    Explicitly excluding the two problematic artifacts instead of including the three (or four?).
    https://github.com/typelevel/cats/tree/v0.7.2#getting-started
    Re "_2.11", see also: https://github.com/sbt/sbt/issues/1518
     */
    exclude("org.typelevel", "cats-laws_2.11")
    exclude("org.typelevel", "cats-kernel-laws_2.11")
    )

  Seq(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-stream" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV,
    "com.typesafe.akka" %% "akka-http" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-testkit" % akkaHttpV,
    "com.iheart" %% "ficus" % "1.4.0",
    "org.webjars" %  "swagger-ui" % "2.1.1",
    "org.scalatest" %% "scalatest" % scalaTestV % Test,
    "io.swagger" % "swagger-parser" % "1.0.22" % Test,
    "org.yaml" % "snakeyaml" % "1.17" % Test,
    "org.broadinstitute" %% "lenthall" % lenthallV
  ) ++ catsDependencies
}

Revolver.settings
resolvers ++= List(
  "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/",
  "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/"
)
