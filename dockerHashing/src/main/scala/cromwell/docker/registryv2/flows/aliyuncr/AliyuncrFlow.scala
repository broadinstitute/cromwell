package cromwell.docker.registryv2.flows.aliyuncr

//import akka.actor.Scheduler
//import akka.http.scaladsl.model.headers.{Authorization, BasicHttpCredentials}
// import akka.http.scaladsl.model._
import akka.stream.ActorMaterializer
//import cromwell.core.DockerCredentials
import AliyunCr._
//import cromwell.docker.DockerHashActor.DockerHashContext
import cromwell.docker.DockerImageIdentifierWithoutHash
//import cromwell.docker.registryv2.flows.aliyuncr.AliyunCrAbstractFlow
import cromwell.docker.registryv2.flows.aliyuncr.AliyunCrAbstractFlow.HttpDockerFlow
import mouse.all._

import scala.concurrent.ExecutionContext

class AliyunCrFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer/*, notUsed: Scheduler*/)
  extends AliyunCrAbstractFlow(httpClientFlow)(ec, materializer/*, notUsed*/) {

  // override val registryHostName = RegistryHostName
  // override val authorizationServerHostName = AuthorizationServerHostName
  // override val serviceName = ServiceName

  // //protected def regionId = 

  // /**
  //   * Builds the list of headers for the token request
  //   */
  //  def buildTokenRequestHeaders(dockerHashContext: DockerHashContext) = {
  //   dockerHashContext.credentials collect {
  //     case DockerCredentials(token, _, _) =>
  //       Authorization(BasicHttpCredentials(token))
  //   }
  // }

  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash): Boolean =
    dockerImageIdentifierWithoutHash.host |> isValidAliyunCrHost
}
