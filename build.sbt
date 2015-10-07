import com.typesafe.sbt.SbtGit._

name := "lenthall"

version := "0.13-SNAPSHOT"

organization := "org.broadinstitute"

scalaVersion := "2.11.7"

libraryDependencies ++= Seq(
  "com.typesafe" % "config" % "1.2.1",
  "org.slf4j" % "slf4j-api" % "1.7.7",
  //---------- Provided libraries -------------------//
  "org.scalaz" %% "scalaz-core" % "7.1.4" % Provided,
  "ch.qos.logback" % "logback-classic" % "1.1.3" % Provided,
  //---------- Test libraries -------------------//
  "org.scalatest" %% "scalatest" % "2.2.5" % Test
)

shellPrompt := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value)}

scalacOptions ++= Seq("-deprecation", "-unchecked", "-feature")

/*
SLF4J initializes itself upon the first logging call.  Because sbt runs tests in parallel it is likely that a second
thread will invoke a second logging call before SLF4J has completed initialization from the first thread's logging call,
leading to these messages:
  SLF4J: The following loggers will not work because they were created
  SLF4J: during the default configuration phase of the underlying logging system.
  SLF4J: See also http://www.slf4j.org/codes.html#substituteLogger
  SLF4J: com.imageworks.common.concurrent.SingleThreadInfiniteLoopRunner

As a workaround, load SLF4J's root logger before starting the unit tests

Source: https://github.com/typesafehub/scalalogging/issues/23#issuecomment-17359537
References:
  http://stackoverflow.com/a/12095245
  http://jira.qos.ch/browse/SLF4J-167
  http://jira.qos.ch/browse/SLF4J-97
*/
testOptions in Test += Tests.Setup(classLoader =>
  classLoader
    .loadClass("org.slf4j.LoggerFactory")
    .getMethod("getLogger", classLoader.loadClass("java.lang.String"))
    .invoke(null, "ROOT")
)
