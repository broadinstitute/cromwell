import com.typesafe.sbt.SbtGit._
import sbt.Keys._
import sbt._

object Version {
  // Upcoming release, or current if we're on a master / hotfix branch
  val cromwellVersion = "24"

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

  val writeVersionConf: Def.Initialize[Task[Seq[File]]] = Def.task {
    val file = (resourceManaged in Compile).value / "cromwell-version.conf"
    val contents =
      s"""|version {
          |  cromwell: "${version.value}"
          |}
          |""".stripMargin
    IO.write(file, contents)
    Seq(file)
  }

  val versionConfCompileSettings = List(resourceGenerators in Compile <+= writeVersionConf)

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

    // For now, obfuscate SNAPSHOTs from sbt's developers: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
    if (isSnapshot) s"$version-SNAP" else version
  }
}
