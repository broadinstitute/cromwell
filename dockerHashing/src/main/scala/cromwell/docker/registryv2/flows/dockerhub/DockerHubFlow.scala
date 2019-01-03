package cromwell.docker.registryv2.flows.dockerhub

import akka.actor.Scheduler
import akka.http.scaladsl.model.headers.{Authorization, BasicHttpCredentials}
import akka.stream.ActorMaterializer
import cromwell.core.DockerCredentials
import DockerHub._
import cromwell.docker.DockerHashActor.DockerHashContext
import cromwell.docker.DockerImageIdentifierWithoutHash
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow
import mouse.all._

import scala.concurrent.ExecutionContext

class DockerHubFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer, scheduler: Scheduler)
  extends DockerRegistryV2AbstractFlow(httpClientFlow)(ec, materializer, scheduler) {

  override val registryHostName = RegistryHostName
  override val authorizationServerHostName = AuthorizationServerHostName
  override val serviceName = ServiceName

  /**
    * Builds the list of headers for the token request
    */
   def buildTokenRequestHeaders(dockerHashContext: DockerHashContext) = {
    dockerHashContext.credentials collect {
      case DockerCredentials(token) =>
        Authorization(BasicHttpCredentials(token))
    }
  }

  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash): Boolean =
    dockerImageIdentifierWithoutHash.host |> isValidDockerHubHost
}
