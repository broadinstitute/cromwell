package cromwell.util.docker

sealed trait DockerIdentifier {
  def name: String

  def registry: DockerRegistry
}

case class DockerDigestIdentifier(name: String, digest: String,
                                  registry: DockerRegistry = DockerRegistry.DockerHub
                                 ) extends DockerIdentifier

case class DockerTagIdentifier(name: String, tag: String = "latest",
                               registry: DockerRegistry = DockerRegistry.DockerHub
                              ) extends DockerIdentifier {
  private lazy val repositoryUri = s"https://${registry.v1Endpoint}/v1/repositories/$name"

  lazy val imagesUri = s"$repositoryUri/images"

  lazy val tagsUri = s"$repositoryUri/tags/$tag"

  lazy val manifestUri = s"https://${registry.v2Endpoint}/v2/$name/manifests/$tag"

  def tokenRequestQueryMap(tokenRequest: DockerV2TokenRequest) = {
    val service = tokenRequest.service
    val scope = tokenRequest.scope getOrElse s"repository:$name:pull"
    Map("service" -> service, "scope" -> scope)
  }
}
