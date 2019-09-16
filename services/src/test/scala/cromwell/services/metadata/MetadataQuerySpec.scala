package cromwell.services.metadata

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.TestProbe
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.services.metadata.MetadataQuerySpec.{CannedResponseReadMetadataWorker, MetadataServiceActor_CustomizeRead}
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataServiceResponse, QueryForWorkflowsMatchingParameters, WorkflowQueryResponse, WorkflowQuerySuccess}
import cromwell.services.metadata.impl.{MetadataServiceActor, MetadataServiceActorSpec}
import org.scalatest.{FlatSpecLike, Matchers}

class MetadataQuerySpec extends TestKitSuite("MetadataQuerySpec") with FlatSpecLike with Matchers  {

  it should "correctly forward requests to read workers and responses back to requesters" in {

    val request = QueryForWorkflowsMatchingParameters(
      parameters = List(("paramName1", "paramValue1"))
    )

    val response = WorkflowQuerySuccess(
      response = WorkflowQueryResponse(Seq.empty, 0),
      meta = None
    )

    val requester = TestProbe("MetadataServiceClientProbe")
    def readWorkerProps() = Props(new CannedResponseReadMetadataWorker(Map(request -> response)))
    val serviceRegistry = TestProbe("ServiceRegistryProbe")
    val metadataService = system.actorOf(MetadataServiceActor_CustomizeRead.props(readWorkerProps, serviceRegistry), "MetadataServiceUnderTest")

    requester.send(metadataService, request)
    requester.expectMsg(response)

  }

}


object MetadataQuerySpec {
  final class MetadataServiceActor_CustomizeRead(config: Config, serviceRegistryActor: ActorRef, readWorkerMaker: () => Props)
    extends MetadataServiceActor(MetadataServiceActorSpec.globalConfigToMetadataServiceConfig(config), config, serviceRegistryActor) {

    override def readMetadataWorkerActorProps(): Props = readWorkerMaker.apply.withDispatcher(cromwell.core.Dispatcher.ServiceDispatcher)
  }

  object MetadataServiceActor_CustomizeRead {
    val config = ConfigFactory.parseString(MetadataServiceActorSpec.ConfigWithoutSummarizer)

    def props(readActorProps: () => Props, serviceRegistryProbe: TestProbe) =
      Props(new MetadataServiceActor_CustomizeRead(config, serviceRegistryProbe.ref, readActorProps))
  }


  final class CannedResponseReadMetadataWorker(cannedResponses: Map[MetadataReadAction, MetadataServiceResponse]) extends Actor {
    override def receive = {
      case msg: MetadataReadAction => sender ! cannedResponses.getOrElse(msg, throw new Exception(s"Unexpected inbound message: $msg"))
    }
  }
}
