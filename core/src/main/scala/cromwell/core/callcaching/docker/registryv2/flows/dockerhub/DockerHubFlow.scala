package cromwell.core.callcaching.docker.registryv2.flows.dockerhub

import akka.http.scaladsl.model.headers.{Authorization, BasicHttpCredentials}
import akka.http.scaladsl.model.{HttpMethod, HttpMethods}
import akka.stream.ActorMaterializer
import cromwell.core.DockerCredentials
import cromwell.core.callcaching.docker.DockerHashActor.DockerHashContext
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow

import scala.concurrent.ExecutionContext

class DockerHubFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer) 
  extends DockerRegistryV2AbstractFlow(httpClientFlow)(ec, materializer) {
  
  override val registryHostName = "registry-1.docker.io"
  override val authorizationServerHostName = "auth.docker.io"
  override val serviceName = Option("registry.docker.io")

  /**
    * Builds the list of headers for the token request
    */
   def buildTokenRequestHeaders(dockerHashContext: DockerHashContext) = {
    dockerHashContext.credentials collect {
      case DockerCredentials(_, token) => Authorization(BasicHttpCredentials(token))
    }
  }

   override def manifestRequestHttpMethod: HttpMethod = HttpMethods.HEAD
}
