import com.typesafe.sbt.SbtGit.GitCommand

name := "lenthall"

organization := "org.broadinstitute"

scalaVersion := "2.11.8"

// Upcoming release, or current if we're on the master branch
git.baseVersion := "0.19"

// Shorten the git commit hash
git.gitHeadCommit := git.gitHeadCommit.value map { _.take(7) }

// Travis will deploy tagged releases, add -SNAPSHOT for all local builds
git.gitUncommittedChanges := true

versionWithGit

val sprayV = "1.3.3"

val akkaV = "2.4.10"

libraryDependencies ++= Seq(
  "com.typesafe" % "config" % "1.3.0",
  "org.slf4j" % "slf4j-api" % "1.7.21",
  //---------- Provided libraries -------------------//
  "org.typelevel" %% "cats" % "0.7.2" % Provided,
  "ch.qos.logback" % "logback-classic" % "1.1.7" % Provided,
  "org.webjars" % "swagger-ui" % "2.2.2" % Provided,
  "io.spray" %% "spray-routing" % sprayV % Provided,
  "io.spray" %% "spray-http" % sprayV % Provided,
  "io.spray" %% "spray-can" % sprayV % Provided,
  "com.typesafe.akka" %% "akka-actor" % akkaV % Provided,
  //---------- Test libraries -------------------//
  "io.spray" %% "spray-testkit" % sprayV % Test,
  "org.scalatest" %% "scalatest" % "3.0.0" % Test
)

shellPrompt := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), git.baseVersion.value)}

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

testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI")
