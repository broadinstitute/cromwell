import com.typesafe.sbt.GitPlugin.autoImport._
import sbt.Keys._

name := "wdl4s"

organization := "org.broadinstitute"

scalaVersion := "2.11.8"

// Upcoming release, or current if we're on the master branch
git.baseVersion := "0.6"

// Shorten the git commit hash
git.gitHeadCommit := git.gitHeadCommit.value map { _.take(7) }

// Travis will deploy tagged releases, add -SNAPSHOT for all local builds
git.gitUncommittedChanges := true

versionWithGit

val sprayJsonV = "1.3.2"

libraryDependencies ++= Seq(
  "com.typesafe.scala-logging" %% "scala-logging" % "3.4.0",
  "io.spray" %% "spray-json" % sprayJsonV,
  "org.scalaz" %% "scalaz-core" % "7.2.5",
  "commons-codec" % "commons-codec" % "1.10",
  "commons-io" % "commons-io" % "2.5",
  "org.apache.commons" % "commons-lang3" % "3.4",
  "com.github.pathikrit" %% "better-files" % "2.16.0",
  //---------- Test libraries -------------------//
  "org.scalatest" %% "scalatest" % "3.0.0" % Test
)

// The reason why -Xmax-classfile-name is set is because this will fail
// to build on Docker otherwise.  The reason why it's 200 is because it
// fails if the value is too close to 256 (even 254 fails).  For more info:
//
// https://github.com/sbt/sbt-assembly/issues/69
// https://github.com/scala/pickling/issues/10
scalacOptions ++= Seq("-deprecation", "-unchecked", "-feature", "-Xmax-classfile-name", "200")

testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI")
