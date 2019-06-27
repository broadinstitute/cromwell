import sbt.Keys._
import sbt._
import Artifactory._

 object Publishing {
  private val buildTimestamp = System.currentTimeMillis() / 1000

   private def artifactoryResolver(isSnapshot: Boolean): Resolver = {
    val repoType = if (isSnapshot) "snapshot" else "release"
    val repoUrl =
      s"${artifactory}libs-$repoType-local;build.timestamp=$buildTimestamp"
    val repoName = "artifactory-publish"
    repoName at repoUrl
  }

   private val artifactoryCredentials: Credentials = {
    val username = sys.env.getOrElse("ARTIFACTORY_USERNAME", "")
    val password = sys.env.getOrElse("ARTIFACTORY_PASSWORD", "")
    Credentials("Artifactory Realm", artifactoryHost, username, password)
  }

   val publishSettings: Seq[Setting[_]] =
  //we only publish to libs-release-local because of a bug in sbt that makes snapshots take
  //priority over the local package cache. see here: https://github.com/sbt/sbt/issues/2687#issuecomment-236586241
    Seq(
      publishTo := Option(artifactoryResolver(false)),
      publishArtifact in Compile := true,
      publishArtifact in Test := true,
      credentials += artifactoryCredentials
    )

   val noPublishSettings: Seq[Setting[_]] =
    Seq(
      publish := {},
      publishLocal := {}
    )
}