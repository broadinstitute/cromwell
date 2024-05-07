import Dependencies._
import com.github.sbt.git.SbtGit._
import sbt.Keys._
import sbt._
import com.github.sbt.git.SbtGit

object Version {
  // Upcoming release, or current if we're on a master / hotfix branch
  val cromwellVersion = "87"

  sealed trait BuildType
  case object Snapshot extends BuildType
  case object Debug extends BuildType
  case object Release extends BuildType
  case object Standard extends BuildType

  def buildType: BuildType = {
    if (isDebug) Debug
    else if (isRelease) Release
    else if (isSnapshot) Snapshot
    else Standard
  }

  /**
    * Returns true if this project should be considered a snapshot.
    *
    * The value is read in directly from the system property `project.isSnapshot` as there were confusing issues with
    * the multi-project and sbt.Keys#isSnapshot().
    *
    * This is the default if no arguments are provided.
    */
  private lazy val isSnapshot: Boolean = getPropOrDefault("project.isSnapshot", default = true)

  /**
    * Returns true if this project should be built in the debugging configuration.
    *
    * Note that this image is much larger than the default build!
    */
  private lazy val isDebug: Boolean = getPropOrDefault("project.isDebug")

  /**
    * Returns `true` if this project should tag a release like `85` in addition to a hash like `85-443a6fc`.
    */
  private lazy val isRelease: Boolean = getPropOrDefault("project.isRelease")

  private def getPropOrDefault(prop: String, default: Boolean = false): Boolean = {
    sys.props.get(prop).map(_.toBoolean).getOrElse(default)
  }

  // Adapted from SbtGit.versionWithGit
  def cromwellVersionWithGit: Seq[Setting[_]] =
    Seq(
      ThisBuild / git.versionProperty := "project.version",
      ThisBuild / git.baseVersion := cromwellVersion,
      ThisBuild / version :=
        makeVersion(versionProperty = git.versionProperty.value,
                    baseVersion = git.baseVersion.?.value,
                    headCommit = git.gitHeadCommit.value
        ),
      ThisBuild / shellPrompt := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), cromwellVersion) }
    )

  val writeProjectVersionConf: Def.Initialize[Task[Seq[File]]] = Def.task {
    writeVersionConf(name.value, (Compile / resourceManaged).value, version.value)
  }

  val writeSwaggerUiVersionConf: Def.Initialize[Task[Seq[File]]] = Def.task {
    writeVersionConf("swagger-ui", (Compile / resourceManaged).value, swaggerUiVersion)
  }

  /**
    * Writes a version.conf compatible with cromwell-common's VersionUtil. Returns the written file wrapped in a Seq to
    * make it compatible for appending to `Compile / resourceGenerators`.
    *
    * Ex:
    * {{{
    * Compile / resourceGenerators += writeVersionConf(name.value, (Compile / resourceManaged).value, version.value)
    * }}}
    *
    * For a project named "my-project", writes a conf named "my-project-version.conf" containing
    * "my.project.version = [version]"
    *
    * @param projectName Name of the project
    * @param directory The managed resource directory, usually `(Compile / resourceManaged).value`
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

  private def makeVersion(versionProperty: String, baseVersion: Option[String], headCommit: Option[String]): String = {
    // The version string passed in via command line settings, if desired.
    def overrideVersion = Option(sys.props(versionProperty))

    val basePrefix = baseVersion.map(_ + "-").getOrElse("")

    // Version string that just uses the commit version.
    def commitVersion: Option[String] = headCommit map (sha => basePrefix + sha.take(7)) // Shorten the git commit hash

    // Version string fallback.
    val unknownVersion = basePrefix + "unknown"

    // Now we fall through the potential version numbers...
    val version = overrideVersion orElse commitVersion getOrElse unknownVersion

    // For now, obfuscate SNAPSHOTs from sbt's developers: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
    // (by calling it `SNAP` instead of `SNAPSHOT`)
    buildType match {
      case Snapshot => s"$version-SNAP"
      case Debug => s"$version-DEBUG"
      case _ => version
    }
  }
}
