name := "Cromwell Backend"

version := "1.0"

scalaVersion := "2.11.7"

organization := ""

scalacOptions := Seq("-unchecked", "-deprecation", "-encoding", "utf8")

resolvers ++= Seq("Typesafe Repository" at "http://repo.typesafe.com/typesafe/releases/")

libraryDependencies ++= {
  val akkaV = "2.3.12"
  Seq("com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % "test",
    "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
    "ch.qos.logback" % "logback-classic" % "1.1.3",
    "org.scalatest" % "scalatest_2.11" % "2.2.5" % "test",
    "com.github.pathikrit" %% "better-files" % "2.13.0",
    "commons-codec" % "commons-codec" % "1.10")
}

assemblyJarName in assembly := "cromwell-backend.jar"