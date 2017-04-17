package cromwell.util.docker

/**
  * A user specification of a docker identifier.
  */
sealed trait DockerIdentifier {
  /**
    * The name of the docker image, including the user/project.
    * Ex:
    * - library/ubuntu
    * - broadinstitute/scala-base-image
    * - broad-dsde-dev/cromwell
    * - broad-dsde-dev/ubuntu
    */
  def name: String

  /** The registry where this image is located. */
  def registry: DockerRegistry
}

/**
  * A user specified image with a digest identifier. NOTE: The digest includes the image type.
  * Ex: broadinstitute/scala-baseimage@sha256:265feb82d1a9fc8593bb1f2605a63cb0e30ad9ac8d1e74d8ed9113bb129c1885
  */
case class DockerDigestIdentifier(name: String, digest: String,
                                  registry: DockerRegistry = DockerRegistry.DockerHub
                                 ) extends DockerIdentifier

/**
  * A user specified image with a tag, defaulting to "latest". The tags may be adjusted to different image versions.
  * Ex:
  * - library/ubuntu:latest
  * - broad-dsde-dev/cromwell:dev
  */
case class DockerTagIdentifier(name: String, tag: String = "latest",
                               registry: DockerRegistry = DockerRegistry.DockerHub
                              ) extends DockerIdentifier {
  private lazy val repositoryUri = s"https://${registry.v1Hostname}/v1/repositories/$name"

  lazy val imagesUri = s"$repositoryUri/images"

  lazy val tagsUri = s"$repositoryUri/tags/$tag"

  lazy val manifestUri = s"https://${registry.v2Hostname}/v2/$name/manifests/$tag"

  def tokenRequestQueryMap(tokenRequest: DockerV2TokenRequest) = {
    val service = tokenRequest.service
    val scope = tokenRequest.scope getOrElse s"repository:$name:pull"
    Map("service" -> service, "scope" -> scope)
  }
}
