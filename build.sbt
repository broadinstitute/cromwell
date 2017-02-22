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
  "org.typelevel" %% "cats" % "0.7.2"
    exclude("org.typelevel", "cats-laws_2.11")
    exclude("org.typelevel", "cats-kernel-laws_2.11"),
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "com.typesafe.akka" %% "akka-http-experimental" % akkaV,
  "com.typesafe.akka" %% "akka-http-spray-json-experimental" % akkaV,
  "com.github.pathikrit" %% "better-files" % "2.13.0",
  //---------- Test libraries -------------------//
  "org.scalatest" %% "scalatest" % "3.0.1" % Test,
  "org.pegdown" % "pegdown" % "1.6.0" % Test
)

testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI", "-h", "target/test-reports")
