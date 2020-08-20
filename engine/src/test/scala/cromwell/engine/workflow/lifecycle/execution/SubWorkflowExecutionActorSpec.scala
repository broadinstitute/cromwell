package cromwell.engine.workflow.lifecycle.execution

import java.util.UUID
import java.util.concurrent.atomic.AtomicInteger

import akka.actor.Props
import akka.testkit.{TestFSMRef, TestProbe}
import com.typesafe.config.ConfigFactory
import cromwell.backend.{AllBackendInitializationData, BackendWorkflowDescriptor, JobExecutionMap}
import cromwell.core._
import cromwell.core.callcaching.CallCachingOff
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder}
import cromwell.database.sql.tables.SubWorkflowStoreEntry
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.execution.SubWorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation.CallPreparationFailed
import cromwell.engine.workflow.lifecycle.execution.job.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import cromwell.engine.workflow.lifecycle.execution.keys.SubWorkflowKey
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.engine.workflow.workflowstore.{RestartableRunning, StartableState, Submitted}
import cromwell.engine.{ContinueWhilePossible, EngineIoFunctions, EngineWorkflowDescriptor}
import cromwell.subworkflowstore.SubWorkflowStoreActor.{QuerySubWorkflow, SubWorkflowFound, SubWorkflowNotFound}
import cromwell.util.WomMocks
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito
import wom.graph.WomIdentifier

import scala.concurrent.duration._
import scala.language.postfixOps

class SubWorkflowExecutionActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Mockito with Eventually {
  
  behavior of "SubWorkflowExecutionActor"

  val serviceRegistryProbe = TestProbe()
  val jobStoreProbe = TestProbe()
  val subWorkflowStoreProbe = TestProbe()
  val callCacheReadActorProbe = TestProbe()
  val callCacheWriteActorProbe = TestProbe()
  val dockerHashActorProbe = TestProbe()
  val ioActorProbe = TestProbe()
  val jobTokenDispenserProbe = TestProbe()
  val preparationActor = TestProbe()
  val subWorkflowActor = TestProbe()
  val deathWatch = TestProbe()
  val parentProbe = TestProbe()
  val parentBackendDescriptor = mock[BackendWorkflowDescriptor]
  val parentWorkflowId: WorkflowId = WorkflowId.randomId()
  parentBackendDescriptor.id returns parentWorkflowId
  val parentWorkflowDescriptor = EngineWorkflowDescriptor(
    WomMocks.mockWorkflowDefinition("workflow"),
    parentBackendDescriptor,
    Map.empty,
    ContinueWhilePossible,
    List.empty,
    CallCachingOff
  )
  val subWorkflow = WomMocks.mockWorkflowDefinition("sub_wf")
  val subWorkflowCall = WomMocks.mockWorkflowCall(WomIdentifier("workflow"), definition = subWorkflow)
  val subKey: SubWorkflowKey = SubWorkflowKey(subWorkflowCall, None, 1)
  val rootConfig = ConfigFactory.load

  val awaitTimeout: FiniteDuration = 10 seconds

  def buildEWEA(startState: StartableState = Submitted) = {
    new TestFSMRef[SubWorkflowExecutionActorState, SubWorkflowExecutionActorData, SubWorkflowExecutionActor](system, Props(
      new SubWorkflowExecutionActor(
        subKey,
        parentWorkflowDescriptor,
        new EngineIoFunctions(List.empty, new AsyncIo(simpleIoActor, DefaultIoCommandBuilder), system.dispatcher),
        Map.empty,
        ioActorProbe.ref,
        serviceRegistryProbe.ref,
        jobStoreProbe.ref,
        subWorkflowStoreProbe.ref,
        callCacheReadActorProbe.ref,
        callCacheWriteActorProbe.ref,
        dockerHashActorProbe.ref,
        jobTokenDispenserProbe.ref,
        BackendSingletonCollection(Map.empty),
        AllBackendInitializationData(Map.empty),
        startState,
        rootConfig,
        new AtomicInteger(),
        fileHashCacheActor = None,
        blacklistCache = None
      ) {
        override def createSubWorkflowPreparationActor(subWorkflowId: WorkflowId) = preparationActor.ref
        override def createSubWorkflowActor(createSubWorkflowActor: EngineWorkflowDescriptor) = subWorkflowActor.ref
      }), parentProbe.ref, s"SubWorkflowExecutionActorSpec-${UUID.randomUUID()}")
  }
  
  it should "Check the sub workflow store when restarting" in {
    val ewea = buildEWEA(startState = RestartableRunning)
    ewea.setState(SubWorkflowPendingState)

    ewea ! Execute
    subWorkflowStoreProbe.expectMsg(QuerySubWorkflow(parentWorkflowId, subKey))
    eventually {
      ewea.stateName shouldBe SubWorkflowCheckingStoreState
    }
  }

  it should "Reuse sub workflow id if found in the store" in {
    import cromwell.core.ExecutionIndex._
    
    val ewea = buildEWEA(startState = RestartableRunning)
    ewea.setState(SubWorkflowCheckingStoreState)
    
    val subWorkflowUuid = WorkflowId.randomId()
    ewea ! SubWorkflowFound(SubWorkflowStoreEntry(Option(0), parentWorkflowId.toString, subKey.node.fullyQualifiedName, subKey.index.fromIndex, subKey.attempt, subWorkflowUuid.toString, None))
    parentProbe.expectMsg(RequestValueStore)
    
    eventually {
      ewea.stateName shouldBe WaitingForValueStore
      ewea.stateData.subWorkflowId shouldBe Some(subWorkflowUuid)
    }
  }

  it should "Fall back to a random Id if the sub workflow id is not found in the store" in {
    val ewea = buildEWEA(startState = RestartableRunning)
    ewea.setState(SubWorkflowCheckingStoreState)

    ewea ! SubWorkflowNotFound(QuerySubWorkflow(parentWorkflowId, subKey))
    parentProbe.expectMsg(RequestValueStore)
    
    eventually {
      ewea.stateName shouldBe WaitingForValueStore
      ewea.stateData.subWorkflowId should not be empty
    }
  }
  
  it should "Request output store" in {
    val ewea = buildEWEA()
    val subWorkflowId = WorkflowId.randomId()
    ewea.setState(WaitingForValueStore, SubWorkflowExecutionActorData(Option(subWorkflowId), None))
    val valueStore = ValueStore.empty
    
    ewea ! valueStore
    preparationActor.expectMsg(CallPreparation.Start(valueStore))
    parentProbe.expectMsg(JobStarting(subKey))
    eventually {
      ewea.stateName shouldBe SubWorkflowPreparingState
    }
  }
  
  it should "Prepare a sub workflow"

  it should "Run a sub workflow" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowPreparingState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId()), None))

    val subWorkflowId = WorkflowId.randomId()
    val subBackendDescriptor = mock[BackendWorkflowDescriptor]
    subBackendDescriptor.id returns subWorkflowId
    val subWorkflowDescriptor = EngineWorkflowDescriptor(
      WomMocks.mockWorkflowDefinition("workflow"),
      subBackendDescriptor,
      Map.empty,
      ContinueWhilePossible,
      List.empty,
      CallCachingOff
    )
    
    ewea ! SubWorkflowPreparationSucceeded(subWorkflowDescriptor, Map.empty)
    subWorkflowActor.expectMsg(WorkflowExecutionActor.ExecuteWorkflowCommand)
    parentProbe.expectMsg(JobRunning(subKey, Map.empty))
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
    val subworkflowId = WorkflowId.randomId()
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(subworkflowId), None))

    deathWatch watch ewea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val outputs: CallOutputs = CallOutputs.empty
    val workflowSuccessfulMessage = WorkflowExecutionSucceededResponse(jobExecutionMap, Set.empty[WorkflowId], outputs)
    ewea ! workflowSuccessfulMessage
    parentProbe.expectMsg(SubWorkflowSucceededResponse(subKey, jobExecutionMap, Set.empty[WorkflowId], outputs))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }

  it should "Relay Workflow Failed message" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId()), None))

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
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId()), None))

    deathWatch watch ewea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val workflowAbortedMessage = WorkflowExecutionAbortedResponse(jobExecutionMap)
    ewea ! workflowAbortedMessage
    parentProbe.expectMsg(SubWorkflowAbortedResponse(subKey, jobExecutionMap))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }

  it should "Relay Workflow Abort command message" in {
    val ewea = buildEWEA()
    ewea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorData(Some(WorkflowId.randomId()), Option(subWorkflowActor.ref)))

    deathWatch watch ewea

    val jobExecutionMap: JobExecutionMap = Map.empty
    ewea ! EngineLifecycleActorAbortCommand
    subWorkflowActor.expectMsg(EngineLifecycleActorAbortCommand)
    subWorkflowActor.reply(WorkflowExecutionAbortedResponse(jobExecutionMap))
    parentProbe.expectMsg(SubWorkflowAbortedResponse(subKey, jobExecutionMap))
    deathWatch.expectTerminated(ewea, awaitTimeout)
  }
}
