package cromwell.engine.workflow.lifecycle.execution

import java.util.UUID

import akka.actor.Props
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend.{AllBackendInitializationData, BackendWorkflowDescriptor, JobExecutionMap}
import cromwell.core._
import cromwell.core.callcaching.CallCachingOff
import cromwell.database.sql.tables.SubWorkflowStoreEntry
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.lifecycle.execution.CallPreparationActor.{CallPreparationFailed, SubWorkflowPreparationSucceeded}
import cromwell.engine.workflow.lifecycle.execution.SubWorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.{ContinueWhilePossible, EngineWorkflowDescriptor}
import cromwell.subworkflowstore.SubWorkflowStoreActor.{QuerySubWorkflow, SubWorkflowFound, SubWorkflowNotFound}
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import wdl4s.{WdlNamespaceWithWorkflow, Workflow, WorkflowCall}

import scala.concurrent.duration._
import scala.language.postfixOps

class SubWorkflowExecutionActorSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito with Eventually {
  
  behavior of "SubWorkflowExecutionActor"

  val serviceRegistryProbe = TestProbe()
  val jobStoreProbe = TestProbe()
  val subWorkflowStoreProbe = TestProbe()
  val callCacheReadActorProbe = TestProbe()
  val jobTokenDispenserProbe = TestProbe()
  val preparationActor = TestProbe()
  val subWorkflowActor = TestProbe()
  val deathWatch = TestProbe()
  val parentProbe = TestProbe()
  val parentBackendDescriptor = mock[BackendWorkflowDescriptor]
  val parentWorkflowId: WorkflowId = WorkflowId.randomId()
  parentBackendDescriptor.id returns parentWorkflowId
  val parentWorkflowDescriptor = EngineWorkflowDescriptor(
    mock[WdlNamespaceWithWorkflow],
    parentBackendDescriptor,
    Map.empty,
    ContinueWhilePossible,
    List.empty,
    CallCachingOff
  )
  val subWorkflow = mock[Workflow]
  subWorkflow.unqualifiedName returns "sub_wf"
  val subWorkflowCall = mock[WorkflowCall]
  subWorkflowCall.fullyQualifiedName returns "foo.bar"
  subWorkflowCall.callable returns subWorkflow
  val subKey = SubWorkflowKey(subWorkflowCall, None, 1)
  
  val awaitTimeout: FiniteDuration = 10 seconds

  def buildEWEA(restart: Boolean = false) = {
    new TestFSMRef[SubWorkflowExecutionActorState, SubWorkflowExecutionActorData, SubWorkflowExecutionActor](system, Props(
      new SubWorkflowExecutionActor(
        subKey,
        WorkflowExecutionActorData.empty(parentWorkflowDescriptor),
        Map.empty,
        serviceRegistryProbe.ref,
        jobStoreProbe.ref,
        subWorkflowStoreProbe.ref,
        callCacheReadActorProbe.ref,
        jobTokenDispenserProbe.ref,
        BackendSingletonCollection(Map.empty),
        AllBackendInitializationData(Map.empty),
        restart
      ) {
        override def createSubWorkflowPreparationActor(subWorkflowId: WorkflowId) = preparationActor.ref
        override def createSubWorkflowActor(createSubWorkflowActor: EngineWorkflowDescriptor) = subWorkflowActor.ref
      }), parentProbe.ref, s"SubWorkflowExecutionActorSpec-${UUID.randomUUID()}")
  }
  
  it should "Check the sub workflow store when restarting" in {
    val ewea = buildEWEA(restart = true)
    ewea.setState(SubWorkflowPendingState)

    ewea ! Execute
    subWorkflowStoreProbe.expectMsg(QuerySubWorkflow(parentWorkflowId, subKey))
    eventually {
      ewea.stateName shouldBe SubWorkflowCheckingStoreState
    }
  }

  it should "Reuse sub workflow id if found in the store" in {
    import cromwell.core.ExecutionIndex._
    
    val ewea = buildEWEA(restart = true)
    ewea.setState(SubWorkflowCheckingStoreState)
    
    val subWorkflowUuid = WorkflowId.randomId()
    ewea ! SubWorkflowFound(SubWorkflowStoreEntry(Option(0), parentWorkflowId.toString, subKey.scope.fullyQualifiedName, subKey.index.fromIndex, subKey.attempt, subWorkflowUuid.toString, None))
    preparationActor.expectMsg(CallPreparationActor.Start)
    parentProbe.expectMsg(JobStarting(subKey))
    
    eventually {
      ewea.stateName shouldBe SubWorkflowPreparingState
      ewea.stateData.subWorkflowId shouldBe Some(subWorkflowUuid)
    }
  }

  it should "Fall back to a random Id if the sub workflow id is not found in the store" in {
    val ewea = buildEWEA(restart = true)
    ewea.setState(SubWorkflowCheckingStoreState)

    ewea ! SubWorkflowNotFound(QuerySubWorkflow(parentWorkflowId, subKey))
    preparationActor.expectMsg(CallPreparationActor.Start)
    parentProbe.expectMsg(JobStarting(subKey))
    
    eventually {
      ewea.stateName shouldBe SubWorkflowPreparingState
      ewea.stateData.subWorkflowId should not be empty
    }
  }
  
  it should "Prepare a sub workflow" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowPendingState)
    
    ewea ! Execute
    preparationActor.expectMsg(CallPreparationActor.Start)
    parentProbe.expectMsg(JobStarting(subKey))
    eventually {
      ewea.stateName shouldBe SubWorkflowPreparingState
    }
  }

  it should "Run a sub workflow" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowPreparingState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId())))

    val subWorkflowId = WorkflowId.randomId()
    val subBackendDescriptor = mock[BackendWorkflowDescriptor]
    subBackendDescriptor.id returns subWorkflowId
    val subWorkflowDescriptor = EngineWorkflowDescriptor(
      mock[WdlNamespaceWithWorkflow],
      subBackendDescriptor,
      Map.empty,
      ContinueWhilePossible,
      List.empty,
      CallCachingOff
    )
    
    ewea ! SubWorkflowPreparationSucceeded(subWorkflowDescriptor, Map.empty)
    subWorkflowActor.expectMsg(WorkflowExecutionActor.ExecuteWorkflowCommand)
    parentProbe.expectMsg(JobRunning(subKey, Map.empty, Option(subWorkflowActor.ref)))
    eventually {
      ewea.stateName shouldBe SubWorkflowRunningState
    }
  }

  it should "Fail a sub workflow if preparation failed" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowPreparingState)
    deathWatch watch ewea

    val subWorkflowKey = mock[SubWorkflowKey]
    val throwable: Exception = new Exception("Expected test exception")
    val preparationFailedMessage: CallPreparationFailed = CallPreparationFailed(subWorkflowKey, throwable)
    ewea ! preparationFailedMessage
    parentProbe.expectMsg(SubWorkflowFailedResponse(subKey, Map.empty, throwable))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }
  
  it should "Relay Workflow Successful message" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId())))

    deathWatch watch ewea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val outputs: CallOutputs = Map.empty[LocallyQualifiedName, JobOutput]
    val workflowSuccessfulMessage = WorkflowExecutionSucceededResponse(jobExecutionMap, outputs)
    ewea ! workflowSuccessfulMessage
    parentProbe.expectMsg(SubWorkflowSucceededResponse(subKey, jobExecutionMap, outputs))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }

  it should "Relay Workflow Failed message" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId())))

    deathWatch watch ewea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val expectedException: Exception = new Exception("Expected test exception")
    
    val workflowSuccessfulMessage = WorkflowExecutionFailedResponse(jobExecutionMap, expectedException)
    ewea ! workflowSuccessfulMessage
    parentProbe.expectMsg(SubWorkflowFailedResponse(subKey, jobExecutionMap, expectedException))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }

  it should "Relay Workflow Aborted message" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId())))

    deathWatch watch ewea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val workflowAbortedMessage = WorkflowExecutionAbortedResponse(jobExecutionMap)
    ewea ! workflowAbortedMessage
    parentProbe.expectMsg(SubWorkflowAbortedResponse(subKey, jobExecutionMap))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }

}
