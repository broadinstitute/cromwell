package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.Done
import akka.actor.CoordinatedShutdown
import com.typesafe.config.ConfigFactory
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.engine.CromwellTerminator
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.WorkflowStoreWriteHeartbeatCommand
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.control.NoStackTrace

class WorkflowStoreHeartbeatWriteActorSpec extends TestKitSuite("WorkflowStoreHeartbeatWriteActorSpec")
  with FlatSpecLike with Matchers with Eventually {

  behavior of "WorkflowStoreHeartbeatWriteActor"

  it should "send a shutdown signal when unable to write workflow heartbeats" in {

    var shutdownCalled = false

    val workflowStore = new InMemoryWorkflowStore {
      override def writeWorkflowHeartbeats(workflowIds: Set[(WorkflowId, OffsetDateTime)],
                                           heartbeatDateTime: OffsetDateTime)
                                          (implicit ec: ExecutionContext): Future[Int] = {
        Future.failed(new RuntimeException("this is expected") with NoStackTrace)
      }
    }

    val workflowStoreAccess = UncoordinatedWorkflowStoreAccess(workflowStore)
    val workflowHeartbeatTypesafeConfig = ConfigFactory.parseString(
      """|danger.debug.only.minimum-heartbeat-ttl = 10 ms
         |system.workflow-heartbeats {
         |  heartbeat-interval = 500 ms
         |  write-batch-size = 1
         |  write-failure-shutdown-duration = 1 s
         |}
         |""".stripMargin
    ).withFallback(ConfigFactory.load())
    val workflowHeartbeatConfig = WorkflowHeartbeatConfig(workflowHeartbeatTypesafeConfig)
    val terminator = new CromwellTerminator {
      override def beginCromwellShutdown(reason: CoordinatedShutdown.Reason): Future[Done] = {
        shutdownCalled = true
        Future.successful(Done)
      }
    }
    val serviceRegistryActor = emptyActor
    val workflowStoreHeartbeatWriteActor = system.actorOf(
      WorkflowStoreHeartbeatWriteActor.props(
        workflowStoreAccess,
        workflowHeartbeatConfig,
        terminator,
        serviceRegistryActor
      ),
      "WorkflowStoreHeartbeatWriteActor"
    )

    val workflowId = WorkflowId.randomId()
    val submissionTime = OffsetDateTime.now()

    implicit val patienceConfig = PatienceConfig(timeout = 30.seconds, interval = 200.milliseconds)
    eventually {
      if (!shutdownCalled) {
        // Workflow heartbeats are not automatically retried. Instead, write the heartbeat each time we retry.
        workflowStoreHeartbeatWriteActor ! WorkflowStoreWriteHeartbeatCommand(workflowId, submissionTime)
        fail("shutdown was not called")
      }
    }
  }
}
