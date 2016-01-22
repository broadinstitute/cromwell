name := "CromwellLogging"

version := "0.1"

scalaVersion := "2.11.7"

libraryDependencies ++= {
  val akkaV = "2.3.12"
  Seq(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % "test",
    "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
    "org.scalatest" % "scalatest_2.11" % "2.2.4" % "test",
    "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
    "ch.qos.logback" % "logback-classic" % "1.1.3",
    // ---- This is only used for Business Logging stuff --//
    "net.logstash.logback" % "logstash-logback-encoder" % "4.5.1"
  )
}


resolvers ++= Seq(
  "Typesafe Repository" at "http://repo.typesafe.com/typesafe/releases/")


