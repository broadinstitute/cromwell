name := "centaur"

version := "1.0"

scalaVersion := "2.11.7"

val sprayV = "1.3.3"
val downgradedSprayV = "1.3.2"
val akkaV = "2.3.14"

resolvers ++= Seq(
   "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/",
   "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/"
)

/***
 * by default log buffering is set to true in sbt, which means
 * that for tests executed in parallel you will not see the 
 * output until the test suite completes.  Setting this to false
 * will not buffer output, but it will be interleaved
 */
// logBuffered in Test := false

libraryDependencies ++= Seq(
  "org.broadinstitute" %% "lenthall" % "0.17",
  "com.typesafe" % "config" % "1.3.0",
  "org.typelevel" %% "cats" % "0.4.1",
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "io.spray" %% "spray-can" % sprayV,
  "io.spray" %% "spray-routing" % sprayV,
  "io.spray" %% "spray-client" % sprayV,
  "io.spray" %% "spray-http" % sprayV,
  "io.spray" %% "spray-json" % downgradedSprayV,
  "org.scalatest" %% "scalatest" % "2.2.6" % Test,
  "com.github.pathikrit" %% "better-files" % "2.13.0"
)
