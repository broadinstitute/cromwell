
import Dependencies._
import Publishing._
import Version._
import sbt.Keys._
import sbt._
import sbtrelease.ReleasePlugin

object Settings {

  val commonResolvers = List(
    Resolver.jcenterRepo,
    "Broad Artifactory Releases" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/"
  )

  name := "wdl4s"

  // The reason why -Xmax-classfile-name is set is because this will fail
  // to build on Docker otherwise.  The reason why it's 200 is because it
  // fails if the value is too close to 256 (even 254 fails).  For more info:
  //
  // https://github.com/sbt/sbt-assembly/issues/69
  // https://github.com/scala/pickling/issues/10
  //
  // Other fancy flags from https://tpolecat.github.io/2017/04/25/scalac-flags.html.
  //
  // Per JG's work in Cromwell, the following can't be turned on without causing piles of errors in wdl4s.  Many of the
  // constructs that are flagged look suspicious and probably warrant further scrutiny, but no time for that now.
  //
  // "-Ywarn-unused:params"              // Warn if a value parameter is unused.

  val baseSettings = List(
    "-unchecked",
    "-deprecation",
    "-feature",
    "-explaintypes",
    "-Xmax-classfile-name", "200",
    "-target:jvm-1.8",
    "-encoding", "UTF-8"
  )

  val warningSettings = List(
    "-Xfuture",
    "-Xlint:adapted-args",
    "-Xlint:by-name-right-associative",
    "-Xlint:constant",
    "-Xlint:delayedinit-select",
    "-Xlint:doc-detached",
    "-Xlint:inaccessible",
    "-Xlint:infer-any",
    "-Xlint:missing-interpolator",
    "-Xlint:nullary-override",
    "-Xlint:nullary-unit",
    "-Xlint:option-implicit",
    "-Xlint:package-object-classes",
    "-Xlint:poly-implicit-overload",
    "-Xlint:private-shadow",
    "-Xlint:stars-align",
    "-Xlint:type-parameter-shadow",
    "-Xlint:unsound-match",
    "-Yno-adapted-args",
    "-Ywarn-dead-code",
    "-Ywarn-numeric-widen",
    "-Ywarn-value-discard",
    "-Ywarn-inaccessible",
    "-Ywarn-unused:implicits",
    "-Ywarn-unused:privates",
    "-Ywarn-unused:locals",
    "-Ywarn-unused:patvars"
  )

  val consoleHostileSettings = List(
    "-Ywarn-unused:imports", // warns about every unused import on every command.
    "-Xfatal-warnings"       // makes those warnings fatal.
  )

  addCompilerPlugin("org.scalamacros" % "paradise" % "2.1.0" cross CrossVersion.full)

  testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI", "-h", "target/test-reports")

  val commonSettings = ReleasePlugin.projectSettings ++ wdl4sVersionWithGit ++ publishingSettings ++ List(
    organization := "org.broadinstitute",
    scalaVersion := "2.12.3",
    resolvers ++= commonResolvers,
    parallelExecution := false,
    scalacOptions ++= (CrossVersion.partialVersion(scalaVersion.value) match {
      case Some((2, 12)) =>
        // The default scalacOptions includes console-hostile options.  These options are overridden specifically below
        // for the `console` target.
        baseSettings ++ warningSettings ++ consoleHostileSettings
      case Some((2, 11)) =>
        // Scala 2.11 takes a simplified set of options
        baseSettings
      case wut => throw new NotImplementedError(s"Found unsupported Scala version $wut. wdl4s does not support versions of Scala other than 2.11 or 2.12.")
    }),
    // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
    scalacOptions in(Compile, doc) := (baseSettings ++ List("-no-link-warnings")),
    // No console-hostile options, otherwise the console is effectively unusable.
    // https://github.com/sbt/sbt/issues/1815
    scalacOptions in(Compile, console) := (baseSettings ++ warningSettings),
    crossScalaVersions := List("2.11.8", "2.12.3")
  )

  val womSettings = List(
    name := "wdl4s-wom",
    libraryDependencies ++= womDependencies
  ) ++ commonSettings

  val wdlSettings = List(
    name := "wdl4s-wdl",
    libraryDependencies ++= wdlDependencies
  ) ++ commonSettings

  val cwlSettings = List(
    name := "wdl4s-cwl",
    libraryDependencies ++= cwlDependencies
  ) ++ commonSettings

  val rootSettings = List(
    name := "wdl4s",
    libraryDependencies ++= rootDependencies
  ) ++ commonSettings

  testOptions in Test += Tests.Argument(TestFrameworks.ScalaTest, "-oDSI", "-h", "target/test-reports")

}
