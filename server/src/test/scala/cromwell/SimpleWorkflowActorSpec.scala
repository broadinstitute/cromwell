
package cromwell

import java.time.OffsetDateTime
import java.util.UUID
import java.util.concurrent.atomic.AtomicInteger

import akka.actor.Props
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.MetadataWatchActor.{FailureMatcher, Matcher}
import cromwell.SimpleWorkflowActorSpec._
import cromwell.core._
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.tokens.DynamicRateLimiter.Rate
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor
import cromwell.engine.workflow.workflowstore.{Submitted, WorkflowHeartbeatConfig, WorkflowToStart}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld.Addressee
import org.scalatest.BeforeAndAfter

import scala.concurrent.duration._
import scala.concurrent.{Await, Promise}

object SimpleWorkflowActorSpec {

  case class TestableWorkflowActorAndMetadataPromise
  (
    workflowActor: TestFSMRef[WorkflowActorState, WorkflowActorData, WorkflowActor],
    supervisor: TestProbe,
    promise: Promise[Unit])
}

class SimpleWorkflowActorSpec extends CromwellTestKitWordSpec with BeforeAndAfter {
  val serviceRegistry = TestProbe().ref

  private def buildWorkflowActor(sampleWdl: SampleWdl,
                                 rawInputsOverride: String,
                                 workflowId: WorkflowId,
                                 matchers: Matcher*): TestableWorkflowActorAndMetadataPromise = {
    val workflowSources = WorkflowSourceFilesWithoutImports(
      workflowSource = Option(sampleWdl.workflowSource()),
      workflowUrl = None,
      workflowRoot = None,
      workflowType = Option("WDL"),
      workflowTypeVersion = None,
      inputsJson = rawInputsOverride,
      workflowOptions = WorkflowOptions.empty,
      labelsJson = "{}",
      warnings = Vector.empty
    )

    val promise = Promise[Unit]()
    val watchActor = system.actorOf(MetadataWatchActor.props(promise, matchers: _*), s"service-registry-$workflowId-${UUID.randomUUID()}")
    val supervisor = TestProbe()
    val config = ConfigFactory.load()
    val workflowToStart = WorkflowToStart(workflowId, OffsetDateTime.now(), workflowSources, Submitted, HogGroup("foo"))
    val callCachingEnabled =false
    val invalidateBadCacheResults =false
    val workflowActor = TestFSMRef(
      factory = new WorkflowActor(workflowToStart, config,
        ioActor = system.actorOf(SimpleIoActor.props),
        callCachingEnabled = callCachingEnabled,
        invalidateBadCacheResults = invalidateBadCacheResults,
        serviceRegistryActor = watchActor,
        workflowLogCopyRouter = system.actorOf(Props.empty, s"workflow-copy-log-router-$workflowId-${UUID.randomUUID()}"),
        jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props),
        subWorkflowStoreActor = system.actorOf(AlwaysHappySubWorkflowStoreActor.props),
        callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props),
        callCacheWriteActor = system.actorOf(EmptyCallCacheWriteActor.props),
        dockerHashActor = system.actorOf(EmptyDockerHashActor.props),
        jobTokenDispenserActor = system.actorOf(JobExecutionTokenDispenserActor.props(serviceRegistry, Rate(100, 1.second), None)),
        backendSingletonCollection = BackendSingletonCollection(Map("Local" -> None)),
        serverMode = true,
        workflowStoreActor = system.actorOf(Props.empty),
        workflowHeartbeatConfig = WorkflowHeartbeatConfig(config),
        totalJobsByRootWf = new AtomicInteger(),
        fileHashCacheActor = None,
        blacklistCache = None),
      supervisor = supervisor.ref,
      name = s"workflow-actor-$workflowId"
    )
    TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise)
  }

  implicit val TestExecutionTimeout = 30.seconds.dilated
  val AwaitAlmostNothing = 100.milliseconds.dilated
  var workflowId: WorkflowId = _
  before {
    workflowId = WorkflowId.randomId()
  }

  "A WorkflowActor" should {
    "start, run, succeed and die" in {
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, _) = buildWorkflowActor(SampleWdl.HelloWorld, SampleWdl.HelloWorld.workflowJson, workflowId)
      val probe = TestProbe()
      probe watch workflowActor
      startingCallsFilter("wf_hello.hello") {
        workflowActor ! StartWorkflowCommand
      }

      probe.expectTerminated(workflowActor, TestExecutionTimeout)
      // Check the parent didn't see anything:
      supervisor.expectNoMessage(AwaitAlmostNothing) // The actor's already terminated. No point hanging around waiting...

    }

    "fail to construct with missing inputs" in {
      val expectedError = "Required workflow input 'wf_hello.hello.addressee' not specified"
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.HelloWorld, "{}", workflowId, failureMatcher)
      val probe = TestProbe()
      probe watch workflowActor
      workflowActor ! StartWorkflowCommand
      Await.result(promise.future, TestExecutionTimeout)
      probe.expectTerminated(workflowActor, AwaitAlmostNothing)
      supervisor.expectMsgPF(AwaitAlmostNothing, "parent should get a failed response") {
        case x: WorkflowFailedResponse =>
          x.workflowId should be(workflowId)
          x.reasons.size should be(1)
          x.reasons.head.getMessage.contains(expectedError) should be(true)
      }
    }

    "fail to construct with inputs of the wrong type" in {
      // TODO WOM: restore offending offensive input name
      val expectedError = "No coercion defined from '3' of type 'spray.json.JsNumber' to 'String'."
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.HelloWorld, s""" { "$Addressee" : 3} """,
        workflowId, failureMatcher)

      val probe = TestProbe()
      probe watch workflowActor
      workflowActor ! StartWorkflowCommand
      try {
        Await.result(promise.future, TestExecutionTimeout)
      } catch {
        case _: Throwable =>
          val info = failureMatcher.nearMissInformation
          fail(s"We didn't see the expected error message $expectedError within $TestExecutionTimeout. ${info.mkString(", ")}")
      }
      probe.expectTerminated(workflowActor, AwaitAlmostNothing)
      supervisor.expectMsgPF(AwaitAlmostNothing, "parent should get a failed response") {
        case x: WorkflowFailedResponse =>
          x.workflowId should be(workflowId)
          x.reasons.size should be(1)
          x.reasons.head.getMessage.contains(expectedError) should be(true)
      }
    }

    "fail when a call fails" in {
      val expectedError = "Job wf_goodbye.goodbye:NA:1 exited with return code 1 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details."
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.GoodbyeWorld, SampleWdl.GoodbyeWorld.workflowJson, workflowId, failureMatcher)
      val probe = TestProbe()
      probe watch workflowActor
      startingCallsFilter("wf_goodbye.goodbye") {
        workflowActor ! StartWorkflowCommand
      }
      Await.result(promise.future, TestExecutionTimeout)
      probe.expectTerminated(workflowActor, 2.seconds)
      supervisor.expectMsgPF(AwaitAlmostNothing, "parent should get a failed response") {
        case x: WorkflowFailedResponse =>
          x.workflowId should be(workflowId)
          x.reasons.size should be(1)
          x.reasons.head.getMessage.contains(expectedError) should be(true)
      }
    }

    "gracefully handle malformed WDL" in {
      val expectedError = "No input bfile found evaluating inputs for expression bfile"
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.CoercionNotDefined, SampleWdl.CoercionNotDefined.workflowJson, workflowId, failureMatcher)
      val probe = TestProbe()
      probe watch workflowActor
      workflowActor ! StartWorkflowCommand
      try {
        Await.result(promise.future, TestExecutionTimeout)
      } catch {
        case _: Throwable =>
          val info = failureMatcher.nearMissInformation
          val errorString =
            if (info.nonEmpty) "We had a near miss: " + info.mkString(", ")
            else s"The expected key was never seen. We saw: [\n  ${failureMatcher.fullEventList.map(e => s"${e.key} -> ${e.value}").mkString("\n  ")}\n]."
          fail(s"We didn't see the expected error message '$expectedError' within $TestExecutionTimeout. $errorString}")
      }
      probe.expectTerminated(workflowActor, AwaitAlmostNothing)
      supervisor.expectMsgPF(AwaitAlmostNothing, "parent should get a failed response") {
        case x: WorkflowFailedResponse =>
          x.workflowId should be(workflowId)
          x.reasons.size should be(1)
          x.reasons.head.getMessage.contains(expectedError) should be(true)
      }
    }
  }

  private def startingCallsFilter[T](callName: String)(block: => T): T = {
    import CromwellTestKitSpec.waitForInfo
    within(TestExecutionTimeout) {
      waitForInfo(s"Starting $callName") {
        block
      }
    }
  }
}

