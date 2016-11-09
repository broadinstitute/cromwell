name := "centaur"

version := "1.0"

scalaVersion := "2.11.8"

val akkaV = "2.4.11" // Note: akka-http branches from akkaV after 2.4.11

/***
 * by default log buffering is set to true in sbt, which means
 * that for tests executed in parallel you will not see the 
 * output until the test suite completes.  Setting this to false
 * will not buffer output, but it will be interleaved
 */
// logBuffered in Test := false

libraryDependencies ++= Seq(
  "com.github.kxbmap" %% "configs" % "0.4.2",
  "com.typesafe" % "config" % "1.3.0",
  "org.typelevel" %% "cats" % "0.7.2",
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "com.typesafe.akka" %% "akka-http-experimental" % akkaV,
  "com.typesafe.akka" %% "akka-http-spray-json-experimental" % akkaV,
  "org.scalatest" %% "scalatest" % "2.2.6" % Test,
  "com.github.pathikrit" %% "better-files" % "2.13.0"
)
