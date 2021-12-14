package cromwell.engine.workflow.lifecycle.execution

import java.util.UUID
import java.util.concurrent.atomic.AtomicInteger

import akka.actor.Props
import akka.testkit.{TestFSMRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.assertion.ManyTimes
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
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataActionAndRespond}
import cromwell.subworkflowstore.SubWorkflowStoreActor.{QuerySubWorkflow, SubWorkflowFound, SubWorkflowNotFound}
import cromwell.util.WomMocks
import org.scalatest.BeforeAndAfterAll
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito
import wom.graph.WomIdentifier

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.control.NoStackTrace

class SubWorkflowExecutionActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Mockito with Eventually with BeforeAndAfterAll {

  behavior of "SubWorkflowExecutionActor"

  var serviceRegistryProbe: TestProbe = _
  var jobStoreProbe: TestProbe = _
  var subWorkflowStoreProbe: TestProbe = _
  var callCacheReadActorProbe: TestProbe = _
  var callCacheWriteActorProbe: TestProbe = _
  var dockerHashActorProbe: TestProbe = _
  var ioActorProbe: TestProbe = _
  var jobRestartCheckTokenDispenserProbe: TestProbe = _
  var jobExecutionTokenDispenserProbe: TestProbe = _
  var preparationActor: TestProbe = _
  var subWorkflowActor: TestProbe = _
  var deathWatch: TestProbe = _
  var parentProbe: TestProbe = _
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

  override def beforeAll(): Unit = {
    serviceRegistryProbe = TestProbe()
    jobStoreProbe = TestProbe()
    subWorkflowStoreProbe = TestProbe()
    callCacheReadActorProbe = TestProbe()
    callCacheWriteActorProbe = TestProbe()
    dockerHashActorProbe = TestProbe()
    ioActorProbe = TestProbe()
    jobExecutionTokenDispenserProbe = TestProbe()
    preparationActor = TestProbe()
    subWorkflowActor = TestProbe()
    deathWatch = TestProbe()
    parentProbe = TestProbe()
  }

  def buildSWEA(startState: StartableState = Submitted) = {
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
        jobRestartCheckTokenDispenserProbe.ref,
        jobExecutionTokenDispenserProbe.ref,
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
    val swea = buildSWEA(startState = RestartableRunning)
    swea.setState(SubWorkflowPendingState)

    swea ! Execute
    subWorkflowStoreProbe.expectMsg(QuerySubWorkflow(parentWorkflowId, subKey))
    eventually {
      swea.stateName shouldBe SubWorkflowCheckingStoreState
    }
  }

  it should "Reuse sub workflow id if found in the store" in {
    import cromwell.core.ExecutionIndex._

    val swea = buildSWEA(startState = RestartableRunning)
    swea.setState(SubWorkflowCheckingStoreState)

    val subWorkflowUuid = WorkflowId.randomId()
    swea ! SubWorkflowFound(SubWorkflowStoreEntry(Option(0), parentWorkflowId.toString, subKey.node.fullyQualifiedName, subKey.index.fromIndex, subKey.attempt, subWorkflowUuid.toString, None))
    parentProbe.expectMsg(RequestValueStore)

    eventually {
      swea.stateName shouldBe WaitingForValueStore
      swea.stateData.subWorkflowId shouldBe Some(subWorkflowUuid)
    }
  }

  it should "Fall back to a random Id if the sub workflow id is not found in the store" in {
    val swea = buildSWEA(startState = RestartableRunning)
    swea.setState(SubWorkflowCheckingStoreState)

    swea ! SubWorkflowNotFound(QuerySubWorkflow(parentWorkflowId, subKey))
    parentProbe.expectMsg(RequestValueStore)

    eventually {
      swea.stateName shouldBe WaitingForValueStore
      swea.stateData.subWorkflowId should not be empty
    }
  }

  it should "Request output store" in {
    val swea = buildSWEA()
    val subWorkflowId = WorkflowId.randomId()
    swea.setState(WaitingForValueStore, SubWorkflowExecutionActorLiveData(Option(subWorkflowId), None))
    val valueStore = ValueStore.empty

    swea ! valueStore
    preparationActor.expectMsg(CallPreparation.Start(valueStore))
    parentProbe.expectMsg(JobStarting(subKey))
    eventually {
      swea.stateName shouldBe SubWorkflowPreparingState
    }
  }

  it should "Prepare a sub workflow"

  it should "Run a sub workflow" in {
    val swea = buildSWEA()
    swea.setState(SubWorkflowPreparingState, SubWorkflowExecutionActorLiveData(Some(WorkflowId.randomId()), None))

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

    swea ! SubWorkflowPreparationSucceeded(subWorkflowDescriptor, Map.empty)
    subWorkflowActor.expectMsg(WorkflowExecutionActor.ExecuteWorkflowCommand)
    parentProbe.expectMsg(JobRunning(subKey, Map.empty))
    eventually {
      swea.stateName shouldBe SubWorkflowRunningState
    }
  }

  it should "Fail a sub workflow if preparation fails" in {
    val swea = buildSWEA()
    swea.setState(SubWorkflowPreparingState, SubWorkflowExecutionActorLiveData(Some(WorkflowId.randomId()), None))
    deathWatch watch swea

    val subWorkflowKey = mock[SubWorkflowKey]
    val throwable: Exception = new Exception("Expected test exception") with NoStackTrace
    val preparationFailedMessage: CallPreparationFailed = CallPreparationFailed(subWorkflowKey, throwable)
    swea ! preparationFailedMessage

    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteSuccess(events)
    }
    parentProbe.expectMsg(SubWorkflowFailedResponse(subKey, Map.empty, throwable))
    deathWatch.expectTerminated(swea, awaitTimeout)
  }

  it should "Relay Workflow Successful message" in {
    val swea = buildSWEA()
    val subworkflowId = WorkflowId.randomId()
    swea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorLiveData(Some(subworkflowId), None))

    deathWatch watch swea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val outputs: CallOutputs = CallOutputs.empty
    val workflowSuccessfulMessage = WorkflowExecutionSucceededResponse(jobExecutionMap, Set.empty[WorkflowId], outputs)
    swea ! workflowSuccessfulMessage
    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteSuccess(events)
    }
    parentProbe.expectMsg(SubWorkflowSucceededResponse(subKey, jobExecutionMap, Set.empty[WorkflowId], outputs))
    deathWatch.expectTerminated(swea, awaitTimeout)
  }

  it should "Relay Workflow Failed message" in {
    val swea = buildSWEA()
    swea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorLiveData(Some(WorkflowId.randomId()), None))

    deathWatch watch swea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val expectedException: Exception = new Exception("Expected test exception") with NoStackTrace

    val workflowFailedMessage = WorkflowExecutionFailedResponse(jobExecutionMap, expectedException)
    swea ! workflowFailedMessage
    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteSuccess(events)
    }
    parentProbe.expectMsg(SubWorkflowFailedResponse(subKey, jobExecutionMap, expectedException))
    deathWatch.expectTerminated(swea, awaitTimeout)
  }

  it should "Switch Succeeded to Failed and try again if the final metadata entry doesn't write" in {
    val swea = buildSWEA()
    swea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorLiveData(Some(WorkflowId.randomId()), None))

    deathWatch watch swea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val expectedException: Exception = new Exception("Expected test exception") with NoStackTrace

    val outputs: CallOutputs = CallOutputs.empty
    val workflowSuccessfulMessage = WorkflowExecutionSucceededResponse(jobExecutionMap, Set.empty[WorkflowId], outputs)
    swea ! workflowSuccessfulMessage
    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteFailure(expectedException, events)
    }

    import ManyTimes.intWithTimes
    10.times {
      serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
        case PutMetadataActionAndRespond(events, _, _) =>
          events.size should be(1)
          events.head.key.key should be("status")
          events.head.value.get.value should be("Failed")
          swea ! MetadataWriteFailure(expectedException, events)
      }
      // Check there are no messages going to the parent yet:
      parentProbe.expectNoMessage(10.millis)
    }

    // Now let's say eventually the write does somehow get through.
    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteSuccess(events)
    }

    // The workflow is now considered failed since lots of metadata is probably lost:
    parentProbe.expectMsgPF(awaitTimeout) {
      case SubWorkflowFailedResponse(`subKey`, `jobExecutionMap`, reason) =>
        reason.getMessage should be("Sub workflow execution actor unable to write final state to metadata")
        reason.getCause should be(expectedException)
    }
    deathWatch.expectTerminated(swea, awaitTimeout)
  }

  it should "Relay Workflow Aborted message" in {
    val swea = buildSWEA()
    swea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorLiveData(Some(WorkflowId.randomId()), None))

    deathWatch watch swea

    val jobExecutionMap: JobExecutionMap = Map.empty
    val workflowAbortedMessage = WorkflowExecutionAbortedResponse(jobExecutionMap)
    swea ! workflowAbortedMessage
    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteSuccess(events)
    }
    parentProbe.expectMsg(SubWorkflowAbortedResponse(subKey, jobExecutionMap))
    deathWatch.expectTerminated(swea, awaitTimeout)
  }

  it should "Relay Workflow Abort command message" in {
    val swea = buildSWEA()
    swea.setState(SubWorkflowRunningState, SubWorkflowExecutionActorLiveData(Some(WorkflowId.randomId()), Option(subWorkflowActor.ref)))

    deathWatch watch swea

    val jobExecutionMap: JobExecutionMap = Map.empty
    swea ! EngineLifecycleActorAbortCommand
    subWorkflowActor.expectMsg(EngineLifecycleActorAbortCommand)
    subWorkflowActor.reply(WorkflowExecutionAbortedResponse(jobExecutionMap))
    serviceRegistryProbe.fishForSpecificMessage(awaitTimeout) {
      case PutMetadataActionAndRespond(events, _, _) => swea ! MetadataWriteSuccess(events)
    }
    parentProbe.expectMsg(SubWorkflowAbortedResponse(subKey, jobExecutionMap))
    deathWatch.expectTerminated(swea, awaitTimeout)
  }
}
