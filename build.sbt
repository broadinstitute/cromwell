import com.typesafe.sbt.GitPlugin.autoImport._
import sbt.Keys._
import sbtassembly.MergeStrategy

name := "wdltool"

organization := "org.broadinstitute"

scalaVersion := "2.12.1"

val wdl4sV = "0.14-7c693a3-SNAP"

lazy val versionSettings = Seq(
  // Upcoming release, or current if we're on the master branch
  git.baseVersion := "0.14",

  // Shorten the git commit hash
  git.gitHeadCommit := git.gitHeadCommit.value map { _.take(7) },

  // Travis will deploy tagged releases, add -SNAPSHOT for all local builds
  git.gitUncommittedChanges := true,

  // For now, obfuscate SNAPSHOTs from sbt's developers: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
  git.uncommittedSignifier := Option("SNAP")
)

versionWithGit ++ versionSettings

assemblyJarName in assembly := "wdltool-" + git.baseVersion.value + ".jar"

logLevel in assembly := Level.Info

resolvers ++= Seq(
  "Broad Artifactory Releases" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
  "Broad Artifactory Snapshots" at "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot/"
)

lazy val catsDependencies = List(
  "org.typelevel" %% "cats" % "0.9.0",
  "org.typelevel" %% "kittens" % "1.0.0-M10",
  "com.github.benhutchison" %% "mouse" % "0.9"
) map (_
  /*
  Exclude test framework cats-laws and its transitive dependency scalacheck.
  If sbt detects scalacheck, it tries to run it.
  Explicitly excluding the two problematic artifacts instead of including the three (or four?).
  https://github.com/typelevel/cats/tree/v0.7.2#getting-started
  Re "_2.12", see also: https://github.com/sbt/sbt/issues/1518
   */
  exclude("org.typelevel", "cats-laws_2.12")
  exclude("org.typelevel", "cats-kernel-laws_2.12")
  )

libraryDependencies ++= Seq(
  "org.broadinstitute" %% "wdl4s" % wdl4sV,
  //---------- Test libraries -------------------//
  "org.scalatest" %% "scalatest" % "3.0.1" % Test
) ++ catsDependencies

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
