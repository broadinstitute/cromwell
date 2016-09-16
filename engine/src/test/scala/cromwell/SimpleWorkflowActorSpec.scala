package cromwell

import java.util.UUID

import akka.actor.Props
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.MetadataWatchActor.{FailureMatcher, Matcher}
import cromwell.SimpleWorkflowActorSpec._
import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld.Addressee
import org.scalatest.BeforeAndAfter

import scala.concurrent.duration._
import scala.concurrent.{Await, Promise}
import scala.language.postfixOps

object SimpleWorkflowActorSpec {

  case class TestableWorkflowActorAndMetadataPromise
  (
    workflowActor: TestFSMRef[WorkflowActorState, WorkflowActorData, WorkflowActor],
    supervisor: TestProbe,
    promise: Promise[Unit])
}

class SimpleWorkflowActorSpec extends CromwellTestkitSpec with BeforeAndAfter {

  private def buildWorkflowActor(sampleWdl: SampleWdl,
                                 rawInputsOverride: String,
                                 workflowId: WorkflowId,
                                 matchers: Matcher*): TestableWorkflowActorAndMetadataPromise = {
    val workflowSources = WorkflowSourceFiles(sampleWdl.wdlSource(), rawInputsOverride, "{}")
    val promise = Promise[Unit]()
    val watchActor = system.actorOf(MetadataWatchActor.props(promise, matchers: _*), s"service-registry-$workflowId-${UUID.randomUUID()}")
    val supervisor = TestProbe()
    val workflowActor = TestFSMRef(
      factory = new WorkflowActor(workflowId, StartNewWorkflow, workflowSources, ConfigFactory.load(),
        serviceRegistryActor = watchActor,
        workflowLogCopyRouter = system.actorOf(Props.empty, s"workflow-copy-log-router-$workflowId-${UUID.randomUUID()}"),
        jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props),
        callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props)),
      supervisor = supervisor.ref,
      name = s"workflow-actor-$workflowId"
    )
    TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise)
  }

  implicit val TestExecutionTimeout = 10.seconds.dilated
  val AwaitAlmostNothing = 100.milliseconds.dilated
  var workflowId: WorkflowId = _
  before {
    workflowId = WorkflowId.randomId()
  }

  "A WorkflowActor" should {
    "start, run, succeed and die" in {
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, _) = buildWorkflowActor(SampleWdl.HelloWorld, SampleWdl.HelloWorld.wdlJson, workflowId)
      val probe = TestProbe()
      probe watch workflowActor
      startingCallsFilter("hello.hello") {
        workflowActor ! StartWorkflowCommand
      }

      probe.expectTerminated(workflowActor, TestExecutionTimeout)
      // Check the parent didn't see anything:
      supervisor.expectNoMsg(AwaitAlmostNothing) // The actor's already terminated. No point hanging around waiting...

    }

    "fail to construct with missing inputs" in {
      val expectedError = "Required workflow input 'hello.hello.addressee' not specified."
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
      val expectedError = "Could not coerce value for 'hello.hello.addressee' into: WdlStringType"
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.HelloWorld, s""" { "$Addressee" : 3} """,
        workflowId, failureMatcher)

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

    "fail when a call fails" in {
      val expectedError = "Call goodbye.goodbye: return code was 1"
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.GoodbyeWorld, SampleWdl.GoodbyeWorld.wdlJson, workflowId, failureMatcher)
      val probe = TestProbe()
      probe watch workflowActor
      startingCallsFilter("goodbye.goodbye") {
        workflowActor ! StartWorkflowCommand
      }
      Await.result(promise.future, TestExecutionTimeout)
      probe.expectTerminated(workflowActor, AwaitAlmostNothing)
      supervisor.expectMsgPF(AwaitAlmostNothing, "parent should get a failed response") {
        case x: WorkflowFailedResponse =>
          x.workflowId should be(workflowId)
          x.reasons.size should be(1)
          x.reasons.head.getMessage.contains(expectedError) should be(true)
      }
    }

    "gracefully handle malformed WDL" in {
      val expectedError = "Input evaluation for Call test1.summary failedVariable 'Can't find bfile' not found"
      val failureMatcher = FailureMatcher(expectedError)
      val TestableWorkflowActorAndMetadataPromise(workflowActor, supervisor, promise) = buildWorkflowActor(SampleWdl.CoercionNotDefined, SampleWdl.CoercionNotDefined.wdlJson, workflowId, failureMatcher)
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
  }

  private def startingCallsFilter[T](callNames: String*)(block: => T): T = {
    import CromwellTestkitSpec.waitForInfo
    within(TestExecutionTimeout) {
      waitForInfo(s"Starting calls: ${callNames.mkString("", ":NA:1, ", ":NA:1")}$$", 1) {
        block
      }
    }
  }
}
