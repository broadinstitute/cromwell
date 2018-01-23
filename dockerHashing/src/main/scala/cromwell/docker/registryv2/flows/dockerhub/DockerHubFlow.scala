package cromwell.docker.registryv2.flows.dockerhub

import akka.actor.Scheduler
import akka.http.scaladsl.model.headers.{Authorization, BasicHttpCredentials}
import akka.stream.ActorMaterializer
import cromwell.core.DockerCredentials
import cromwell.docker.DockerHashActor.DockerHashContext
import cromwell.docker.DockerImageIdentifierWithoutHash
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow
import cromwell.docker.registryv2.flows.dockerhub.DockerHubFlow._

import scala.concurrent.ExecutionContext

object DockerHubFlow {
  private val DockerHubDefaultRegistry = "index.docker.io"
  private val RegistryHostName = "registry-1.docker.io"
  private val AuthorizationServerHostName = "auth.docker.io"
  private val DockerHubImplicitRegistryHostName = "docker.io"
  private val ServiceName = Option("registry.docker.io")
  private val validDockerHubHosts = Array(DockerHubDefaultRegistry, RegistryHostName, DockerHubImplicitRegistryHostName)
}

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
      case DockerCredentials(_, token) => 
        Authorization(BasicHttpCredentials(token))
    }
  }
  
  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash) = {
    dockerImageIdentifierWithoutHash.host match {
        // If no host is provided, assume dockerhub
      case None => true
      case Some(host) if validDockerHubHosts.contains(host) => true
      case _ => false
    }
  }
}
