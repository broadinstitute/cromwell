import com.typesafe.sbt.SbtGit._
import sbt.Keys._
import sbt._

object Version {
  // Upcoming release, or current if we're on the master branch
  val cromwellVersion = "0.21"

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

    // The project isSnapshot string passed in via command line settings, if desired.
    val isSnapshot = sys.props.get("project.isSnapshot").forall(_.toBoolean)

    if (isSnapshot) s"$version-SNAPSHOT" else version
  }
}
