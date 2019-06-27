import scala.sys.process._

 object Version {

   def createVersion(baseVersion: String) = {
    def getLastCommitFromGit = { s"""git rev-parse --short HEAD""" !! }

     // either specify git hash as an env var or derive it
    // if building from the broadinstitute/scala-baseimage docker image use env var
    // (scala-baseimage doesn't have git in it)
    val lastCommit = sys.env.getOrElse("GIT_HASH", getLastCommitFromGit ).trim()
    val version = baseVersion + "-" + lastCommit

     // The project isSnapshot string passed in via command line settings, if desired.
    val isSnapshot = sys.props.getOrElse("project.isSnapshot", "true").toBoolean

     // For now, obfuscate SNAPSHOTs from sbt's developers: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
    if (isSnapshot) s"$version-SNAP" else version
  }
}