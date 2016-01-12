name := "Cromwell Backend"

version := "1.0"

scalaVersion := "2.11.7"

organization := ""

scalacOptions := Seq("-unchecked", "-deprecation", "-encoding", "utf8")

resolvers ++= Seq("Typesafe Repository" at "http://repo.typesafe.com/typesafe/releases/")

libraryDependencies ++= {
  val akkaV = "2.3.12"
  Seq("com.typesafe.akka" %% "akka-actor" % akkaV)
}

assemblyJarName in assembly := "cromwell-backend.jar"
