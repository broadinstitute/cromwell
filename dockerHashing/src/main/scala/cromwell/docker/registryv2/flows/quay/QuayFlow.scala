package cromwell.docker.registryv2.flows.quay

import akka.actor.Scheduler
import akka.http.scaladsl.model.HttpHeader
import akka.stream.scaladsl.{Flow, GraphDSL, Source}
import akka.stream.{ActorMaterializer, FanOutShape2}
import cromwell.docker.DockerHashActor.{DockerHashContext, DockerHashResponse}
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow

import scala.concurrent.ExecutionContext

class QuayFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer, scheduler: Scheduler) extends DockerRegistryV2AbstractFlow(httpClientFlow) {
  override protected def registryHostName: String = "quay.io"
  // Not used for now because we bypass the token part as quay doesn't require one for public images 
  override protected def authorizationServerHostName: String = "quay.io"
  // Not used for now, same reason as above
  override protected def buildTokenRequestHeaders(dockerHashContext: DockerHashContext): List[HttpHeader] = List.empty

  override private [registryv2] val tokenFlow = GraphDSL.create() { implicit b =>
    // Flow that always returns a successful empty token
    val noTokenFlow = b.add(Flow[DockerHashContext].map((None, _)))
    // this token flow never fails as it returns an empty token. 
    // However it still needs a failure port. This creates an empty source that will act as the failure output port.
    val emptyFailurePort = b.add(Source.empty[(DockerHashResponse, DockerHashContext)])
    new FanOutShape2(noTokenFlow.in, noTokenFlow.out, emptyFailurePort.out)
  }
}
