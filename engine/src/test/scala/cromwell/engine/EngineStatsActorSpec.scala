package cromwell.engine

import akka.testkit.{ImplicitSender, TestActorRef}
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core._
import cromwell.engine.EngineStatsActor.{EngineStats, StatsQuery}
import cromwell.engine.workflow.WorkflowActor.WorkflowStateTransition
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionAbortedResponse, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionDiff
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.duration._

class EngineStatsActorSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito with ImplicitSender {
  behavior of "EngineStatsActor"
  
  val statsActor = TestActorRef(new EngineStatsActor(Duration.Zero))
  
  val workflowId1 = WorkflowId.randomId()
  val workflowId2 = WorkflowId.randomId()
  
  val jobKey1 = mock[BackendJobDescriptorKey]
  val jobKey2 = mock[BackendJobDescriptorKey]
  val jobKey3 = mock[BackendJobDescriptorKey]
  val jobKey4 = mock[BackendJobDescriptorKey]
  
  it should "update count for new workflows and jobs" in {
    statsActor ! WorkflowExecutionDiff(Map(
      jobKey1 -> ExecutionStatus.NotStarted,
      jobKey2 -> ExecutionStatus.Running,
      jobKey3 -> ExecutionStatus.Running
    ))
    statsActor ! WorkflowStateTransition(workflowId1, WorkflowSubmitted)
    statsActor ! WorkflowStateTransition(workflowId2, WorkflowSubmitted)
    statsActor ! StatsQuery
    
    expectMsg(EngineStats(Map("Submitted" -> 2), Map("Running" -> 2, "NotStarted" -> 1)))
  }

  it should "update count for existing workflows and jobs" in {
    statsActor ! WorkflowExecutionDiff(Map(
      jobKey1 -> ExecutionStatus.Running,
      jobKey2 -> ExecutionStatus.Failed,
      jobKey3 -> ExecutionStatus.Done
    ))
    statsActor ! WorkflowStateTransition(workflowId1, WorkflowRunning)
    statsActor ! StatsQuery

    expectMsg(EngineStats(Map("Submitted" -> 1, "Running" -> 1), Map("Running" -> 1, "Failed" -> 1, "Done" -> 1)))
  }

  it should "remove workflows from count when workflow completes" in {
    statsActor ! WorkflowStateTransition(workflowId1, WorkflowSucceeded)
    statsActor ! StatsQuery
    expectMsg(EngineStats(Map("Submitted" -> 1), Map("Running" -> 1, "Failed" -> 1, "Done" -> 1)))
    
    statsActor ! WorkflowStateTransition(workflowId2, WorkflowFailed)
    statsActor ! StatsQuery
    expectMsg(EngineStats(Map.empty, Map("Running" -> 1, "Failed" -> 1, "Done" -> 1)))
  }

  it should "remove jobs from count when workflow execution completes" in {
    val weSuccess = WorkflowExecutionSucceededResponse(Map(mock[BackendWorkflowDescriptor] -> List(jobKey1)), Map.empty)
    val weFailure = WorkflowExecutionFailedResponse(Map(mock[BackendWorkflowDescriptor] -> List(jobKey2)), new Exception(""))
    val weAbort = WorkflowExecutionAbortedResponse(Map(mock[BackendWorkflowDescriptor] -> List(jobKey3)))
    
    statsActor ! weSuccess
    statsActor ! StatsQuery
    expectMsg(EngineStats(Map.empty, Map("Failed" -> 1, "Done" -> 1)))
    
    statsActor ! weFailure
    statsActor ! StatsQuery
    expectMsg(EngineStats(Map.empty, Map("Done" -> 1)))
    
    statsActor ! weAbort
    statsActor ! StatsQuery
    expectMsg(EngineStats(Map.empty, Map.empty))
  }
}
