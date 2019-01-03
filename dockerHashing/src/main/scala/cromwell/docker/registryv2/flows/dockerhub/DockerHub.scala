package cromwell.docker.registryv2.flows.dockerhub

object DockerHub {
  val DockerHubDefaultRegistry = "index.docker.io"
  val RegistryHostName = "registry-1.docker.io"
  val AuthorizationServerHostName = "auth.docker.io"
  val DockerHubImplicitRegistryHostName = "docker.io"
  val ServiceName = Option("registry.docker.io")
  val validDockerHubHosts = Array(DockerHubDefaultRegistry, RegistryHostName, DockerHubImplicitRegistryHostName)

  def isValidDockerHubHost(host: Option[String]): Boolean =
    host match {
      // If no host is provided, assume dockerhub
      case None => true
      case Some(h) if validDockerHubHosts.contains(h) => true
      case _ => false
    }
}
