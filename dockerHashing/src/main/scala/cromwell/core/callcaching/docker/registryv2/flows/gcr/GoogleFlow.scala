package cromwell.core.callcaching.docker.registryv2.flows.gcr

import akka.actor.Scheduler
import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition}
import akka.stream.{ActorMaterializer, FlowShape, ThrottleMode}
import cromwell.core.callcaching.docker.DockerHashActor.{DockerHashContext, DockerHashResponse, DockerHashUnknownRegistry}
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow
import cromwell.core.callcaching.docker.{DockerFlow, DockerImageIdentifierWithoutHash}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

class GoogleFlow(httpClientFlow: HttpDockerFlow, queriesPer100Sec: Int)(implicit ec: ExecutionContext, materializer: ActorMaterializer, scheduler: Scheduler) extends DockerFlow {
  
  private val gcrFlow = new GcrFlow(httpClientFlow)
  private val usGcrFlow = new GcrUsFlow(httpClientFlow)
  private val euGcrFlow = new GcrEuFlow(httpClientFlow)
  
  private val googleFlows = List(gcrFlow, usGcrFlow, euGcrFlow)
  
  // Throttle all requests to GCR
  private val throttler = Flow[DockerHashContext].throttle(queriesPer100Sec, 100 seconds, queriesPer100Sec, ThrottleMode.Shaping)

  override def buildFlow() = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // This flow should only receive images that have been vetted by accepts.
    def mapContextToFlow(dockerHashContext: DockerHashContext) = {
      googleFlows.indexWhere(_.accepts(dockerHashContext.dockerImageID)) match {
        case -1 => googleFlows.length
        case found => found
      }
    }
    
    // If we ever get an image that none of the above flow can process,
    // which shouldn't happen because the caller should have called "accepts" before
    // to make sure this flow can handle the image
    val failFlow = builder.add(Flow[DockerHashContext].map(dockerContext => (DockerHashUnknownRegistry(dockerContext.request), dockerContext)))
    
    val partition = builder.add(Partition[DockerHashContext](googleFlows.length + 1, mapContextToFlow))
    val merge = builder.add(Merge[(DockerHashResponse, DockerHashContext)](googleFlows.length + 1))
    val throttle = builder.add(throttler)
    
    throttle ~> partition
                googleFlows.indices foreach { i => partition.out(i) ~> googleFlows(i).buildFlow() ~> merge.in(i) }
                partition.out(googleFlows.length) ~> failFlow ~> merge.in(googleFlows.length)
    
    FlowShape(throttle.in, merge.out)
  }
  
  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash) = {
    googleFlows.exists(_.accepts(dockerImageIdentifierWithoutHash))
  }
}
