name := "centaur"

version := "1.0"

scalaVersion := "2.12.1"

val akkaV = "2.4.17"
val akkaHttpV = "10.0.5"

/***
 * by default log buffering is set to true in sbt, which means
 * that for tests executed in parallel you will not see the 
 * output until the test suite completes.  Setting this to false
 * will not buffer output, but it will be interleaved
 */
// logBuffered in Test := false

libraryDependencies ++= Seq(
  "com.github.kxbmap" %% "configs" % "0.4.4",
  "com.typesafe" % "config" % "1.3.0",
  "org.typelevel" %% "cats" % "0.9.0"
    exclude("org.typelevel", "cats-laws_2.11")
    exclude("org.typelevel", "cats-kernel-laws_2.11"),
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "com.typesafe.akka" %% "akka-http" % akkaHttpV,
  "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
  "com.github.pathikrit" %% "better-files" % "3.0.0",
  //---------- Test libraries -------------------//
  "org.scalatest" %% "scalatest" % "3.0.1" % Test,
  "org.pegdown" % "pegdown" % "1.6.0" % Test
)

testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI", "-h", "target/test-reports")
