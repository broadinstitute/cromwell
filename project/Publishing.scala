import sbt.Keys._
import sbt._

object Publishing {
  private val buildTimestamp = System.currentTimeMillis() / 1000

  private def artifactoryResolver(isSnapshot: Boolean): Resolver = {
    val repoType = if (isSnapshot) "snapshot" else "release"
    val repoUrl =
      s"https://artifactory.broadinstitute.org/artifactory/libs-$repoType-local;build.timestamp=$buildTimestamp"
    val repoName = "artifactory-publish"
    repoName at repoUrl
  }

  private val artifactoryCredentials: Credentials = {
    val username = sys.env.getOrElse("ARTIFACTORY_USERNAME", "")
    val password = sys.env.getOrElse("ARTIFACTORY_PASSWORD", "")
    Credentials("Artifactory Realm", "artifactory.broadinstitute.org", username, password)
  }

  def publishingSettings: Seq[Setting[_]] =
    Seq(publishTo := Option(artifactoryResolver(isSnapshot.value)), credentials += artifactoryCredentials)
}
