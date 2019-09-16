import Dependencies._
import com.typesafe.sbt.SbtGit._
import sbt.Keys._
import sbt._

object Version {
  // Upcoming release, or current if we're on a master / hotfix branch
  val cromwellVersion = "46"

  /**
    * Returns true if this project should be considered a snapshot.
    *
    * The value is read in directly from the system property `project.isSnapshot` as there were confusing issues with
    * the multi-project and sbt.Keys#isSnapshot().
    */
  val isSnapshot = sys.props.get("project.isSnapshot").forall(_.toBoolean)

  // Adapted from SbtGit.versionWithGit
  def cromwellVersionWithGit: Seq[Setting[_]] =
    Seq(
      git.versionProperty in ThisBuild := "project.version",
      git.baseVersion in ThisBuild := cromwellVersion,
      version in ThisBuild :=
        makeVersion(
          versionProperty = git.versionProperty.value,
          baseVersion = git.baseVersion.?.value,
          headCommit = git.gitHeadCommit.value),
      shellPrompt in ThisBuild := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), cromwellVersion) }
    )

  val writeProjectVersionConf: Def.Initialize[Task[Seq[File]]] = Def.task {
    writeVersionConf(name.value, (resourceManaged in Compile).value, version.value)
  }

  val writeSwaggerUiVersionConf: Def.Initialize[Task[Seq[File]]] = Def.task {
    writeVersionConf("swagger-ui", (resourceManaged in Compile).value, swaggerUiVersion)
  }

  /**
    * Writes a version.conf compatible with cromwell-common's VersionUtil. Returns the written file wrapped in a Seq to
    * make it compatible for appending to `resourceGenerators in Compile`.
    *
    * Ex:
    * {{{
    * resourceGenerators in Compile += writeVersionConf(name.value, (resourceManaged in Compile).value, version.value)
    * }}}
    *
    * For a project named "my-project", writes a conf named "my-project-version.conf" containing
    * "my.project.version = [version]"
    *
    * @param projectName Name of the project
    * @param directory The managed resource directory, usually `(resourceManaged in Compile).value`
    * @param version The version to write
    * @return The written file
    */
  private def writeVersionConf(projectName: String, directory: File, version: String): Seq[File] = {
    val file = directory / s"$projectName-version.conf"
    val contents =
      s"""|${projectName.replace("-", ".")}.version: "$version"
          |""".stripMargin
    IO.write(file, contents)
    List(file)
  }

  private def makeVersion(versionProperty: String,
                          baseVersion: Option[String],
                          headCommit: Option[String]): String = {
    // The version string passed in via command line settings, if desired.
    def overrideVersion = Option(sys.props(versionProperty))

    val basePrefix = baseVersion.map(_ + "-").getOrElse("")

    // Version string that just uses the commit version.
    def commitVersion: Option[String] = headCommit map (sha => basePrefix + sha.take(7)) // Shorten the git commit hash

    // Version string fallback.
    val unknownVersion = basePrefix + "unknown"

    //Now we fall through the potential version numbers...
    val version = overrideVersion orElse commitVersion getOrElse unknownVersion

    // For now, obfuscate SNAPSHOTs from sbt's developers: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
    if (isSnapshot) s"$version-SNAP" else version
  }
}
