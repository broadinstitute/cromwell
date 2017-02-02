package cromwell.core.callcaching.docker.registryv2.flows.gcr

import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition}
import akka.stream.{ActorMaterializer, FlowShape, ThrottleMode}
import cromwell.core.callcaching.docker.DockerHashActor.{DockerHashContext, DockerHashResponse, DockerHashUnknownRegistry}
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow
import cromwell.core.callcaching.docker.{DockerFlow, DockerImageIdentifierWithoutHash}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

class GoogleFlow(httpClientFlow: HttpDockerFlow, queriesPer100Sec: Int)(implicit ec: ExecutionContext, materializer: ActorMaterializer) extends DockerFlow {
  
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
      if (gcrFlow.accepts(dockerHashContext.dockerImageID)) 0
      else if (usGcrFlow.accepts(dockerHashContext.dockerImageID)) 1
      else if (usGcrFlow.accepts(dockerHashContext.dockerImageID)) 2
      else 3
    }
    
    val gcr = builder.add(gcrFlow.buildFlow())
    val usGcr = builder.add(usGcrFlow.buildFlow())
    val euGcr = builder.add(euGcrFlow.buildFlow())
    // If we ever get an image that none of the above flow can process,
    // which shouldn't happen because the caller should have called "accepts" before
    // to make sure this flow can handle the image
    val failFlow = builder.add(Flow[DockerHashContext].map(dockerContext => (DockerHashUnknownRegistry(dockerContext.dockerImageID), dockerContext)))
    
    val partition = builder.add(Partition[DockerHashContext](4, mapContextToFlow))
    val merge = builder.add(Merge[(DockerHashResponse, DockerHashContext)](4))
    val throttle = builder.add(throttler)
    
    throttle ~> partition
                partition.out(0) ~> gcr ~> merge.in(0)
                partition.out(1) ~> usGcr ~> merge.in(1)
                partition.out(2) ~> euGcr ~> merge.in(2)
                partition.out(3) ~> failFlow ~> merge.in(3)
    
    FlowShape(throttle.in, merge.out)
  }
  
  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash) = {
    googleFlows.exists(_.accepts(dockerImageIdentifierWithoutHash))
  }
}
