package cromwell.docker.registryv2.flows.dockerhub

import cromwell.core.DockerCredentials
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.registryv2.flows.dockerhub.DockerHub._
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import mouse.all._
import org.http4s.headers.Authorization

class DockerHubRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {

  override def registryHostName(dockerImageIdentifier: DockerImageIdentifier) = RegistryHostName
  override def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier) = AuthorizationServerHostName
  override val serviceName = ServiceName

  /**
    * Builds the list of headers for the token request
    */
   override def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext) = {
    dockerInfoContext.credentials collect {
      case DockerCredentials(token) =>
        Authorization(org.http4s.BasicCredentials(token))
    }
  }

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean =
    dockerImageIdentifier.host |> isValidDockerHubHost
}
