import com.typesafe.sbt.GitPlugin.autoImport._
import sbt.Keys._
import sbtassembly.MergeStrategy
import com.typesafe.sbt.SbtGit.GitCommand

name := "wdl4s"

organization := "org.broadinstitute"

scalaVersion := "2.11.7"

// Upcoming release, or current if we're on the master branch
git.baseVersion := "0.3"

// Shorten the git commit hash
git.gitHeadCommit := git.gitHeadCommit.value map { _.take(7) }

// Travis will deploy tagged releases, add -SNAPSHOT for all local builds
git.gitUncommittedChanges := true

versionWithGit

assemblyJarName in assembly := "wdl4s-" + git.baseVersion.value + ".jar"

logLevel in assembly := Level.Info

val DowngradedSprayV = "1.3.1"

libraryDependencies ++= Seq(
  "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
  "io.spray" %% "spray-json" % DowngradedSprayV,
  "org.scalaz" %% "scalaz-core" % "7.1.3",
  "commons-codec" % "commons-codec" % "1.10",
  "commons-io" % "commons-io" % "2.4",
  "org.apache.commons" % "commons-lang3" % "3.4",
  "com.github.pathikrit" %% "better-files" % "2.13.0",
  //---------- Test libraries -------------------//
  "org.scalatest" %% "scalatest" % "2.2.5" % Test
)

val customMergeStrategy: String => MergeStrategy = {
  case x if Assembly.isConfigFile(x) =>
    MergeStrategy.concat
  case PathList(ps@_*) if Assembly.isReadme(ps.last) || Assembly.isLicenseFile(ps.last) =>
    MergeStrategy.rename
  case PathList("META-INF", path@_*) =>
    path map {
      _.toLowerCase
    } match {
      case ("manifest.mf" :: Nil) | ("index.list" :: Nil) | ("dependencies" :: Nil) =>
        MergeStrategy.discard
      case ps@(x :: xs) if ps.last.endsWith(".sf") || ps.last.endsWith(".dsa") =>
        MergeStrategy.discard
      case "plexus" :: xs =>
        MergeStrategy.discard
      case "spring.tooling" :: xs =>
        MergeStrategy.discard
      case "services" :: xs =>
        MergeStrategy.filterDistinctLines
      case ("spring.schemas" :: Nil) | ("spring.handlers" :: Nil) =>
        MergeStrategy.filterDistinctLines
      case _ => MergeStrategy.deduplicate
    }
  case "asm-license.txt" | "overview.html" | "cobertura.properties" =>
    MergeStrategy.discard
  case _ => MergeStrategy.deduplicate
}

assemblyMergeStrategy in assembly := customMergeStrategy

// The reason why -Xmax-classfile-name is set is because this will fail
// to build on Docker otherwise.  The reason why it's 200 is because it
// fails if the value is too close to 256 (even 254 fails).  For more info:
//
// https://github.com/sbt/sbt-assembly/issues/69
// https://github.com/scala/pickling/issues/10
scalacOptions ++= Seq("-deprecation", "-unchecked", "-feature", "-Xmax-classfile-name", "200")
