package cromwell.services.metadata.impl

import java.util.UUID

import akka.actor.Props
import akka.testkit.{EventFilter, TestProbe}
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataService.{DeleteMetadataAction, DeleteMetadataFailedResponse}
import org.scalatest.FlatSpecLike

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._

class DeleteMetadataActorSpec extends TestKitSuite with FlatSpecLike {

  private val workflowId = UUID.randomUUID().toString
  private val deleteMetadataActor = system.actorOf(Props(new DeleteMetadataActor {
    override def deleteNonLabelMetadataEntriesForWorkflow(rootWorkflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Int] = {
      Future.failed(new RuntimeException("Test exception"))
    }
  }))

  it should "try to repeat delete action required amount of times in case of failure" in {
    val probe = TestProbe()
    EventFilter.error(start = "Cannot delete metadata. Remaining number of attempts", occurrences = 4) intercept {
      probe.send(deleteMetadataActor, DeleteMetadataAction(WorkflowId.fromString(workflowId), probe.ref, maxAttempts = 5))
    }
    probe.expectMsgPF(5.seconds) {
      case DeleteMetadataFailedResponse(_, _) => // success
    }
  }
}


