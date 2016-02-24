name := "centaur"

version := "1.0"

scalaVersion := "2.11.7"

val sprayV = "1.3.3"
val downgradedSprayV = "1.3.2"
val akkaV = "2.3.14"

libraryDependencies ++= Seq(
  "com.typesafe" % "config" % "1.3.0",
  "com.typesafe.akka" %% "akka-actor" % akkaV,
  "io.spray" %% "spray-can" % sprayV,
  "io.spray" %% "spray-routing" % sprayV,
  "io.spray" %% "spray-client" % sprayV,
  "io.spray" %% "spray-http" % sprayV,
  "io.spray" %% "spray-json" % downgradedSprayV,
  "org.scalatest" %% "scalatest" % "2.2.6" % Test
)

