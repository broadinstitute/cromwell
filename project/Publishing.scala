import sbt.Keys._
import sbt._

object Publishing {
  private val buildTimestamp = System.currentTimeMillis() / 1000

  private def artifactoryResolver(isSnapshot: Boolean): Resolver = {
    val repoType = if (isSnapshot) "snapshot" else "release"
    val repoUrl =
      s"https://broadinstitute.jfrog.io/broadinstitute/libs-$repoType-local;build.timestamp=$buildTimestamp"
    val repoName = "artifactory-publish"
    repoName at repoUrl
  }

  private val artifactoryCredentials: Seq[Credentials] = {
    val credentialsFile = file("target/ci/resources/artifactory_credentials.properties").getAbsoluteFile
    if (credentialsFile.exists)
      List(Credentials(credentialsFile))
    else
      Nil
  }

  def publishingSettings: Seq[Setting[_]] =
    Seq(publishTo := Option(artifactoryResolver(isSnapshot.value)), credentials ++= artifactoryCredentials)
}
