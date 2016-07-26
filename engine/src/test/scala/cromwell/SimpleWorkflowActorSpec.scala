package cromwell

import java.util.UUID

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.SimpleWorkflowActorSpec._
import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor.{StartNewWorkflow, StartWorkflowCommand}
import cromwell.services.metadata.MetadataService
import MetadataService.PutMetadataAction
import cromwell.jobstore.{JobStoreActor, WriteCountingJobStoreDatabase}
import cromwell.services.metadata.MetadataEvent
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld.Addressee

import scala.concurrent.duration._
import scala.concurrent.{Await, Promise}
import scala.language.postfixOps


object SimpleWorkflowActorSpec {

  val failurePatternMatcher = """failures\[\d*\].message""".r

  object MetadataWatchActor {
    def props(matcher: Option[FailureMatcher], promise: Promise[Unit]): Props = Props(MetadataWatchActor(matcher, promise))
  }

  // This actor stands in for the service registry and watches for metadata messages that match the optional `matcher`.
  // If there is no predicate then use `emptyBehavior` which ignores all messages.  This is here because there is no
  // WorkflowManagerActor in this test that would spin up a real ServiceRegistry.
  case class MetadataWatchActor(matcher: Option[FailureMatcher], promise: Promise[Unit]) extends Actor {
    override def receive = matcher match {
      case None =>
        promise.trySuccess(())
        Actor.emptyBehavior
      case Some(m) => {
        case PutMetadataAction(events) if m.matches(events) =>
          promise.trySuccess(())
          context.stop(self)
      }
    }
  }

  case class FailureMatcher(value: String) {
    def matches(events: Traversable[MetadataEvent]): Boolean = {
      events.exists(e => failurePatternMatcher.findFirstIn(e.key.key).isDefined && e.value.exists { v => v.value.contains(value) })
    }
  }

  case class WorkflowActorAndMetadataPromise(workflowActor: ActorRef, promise: Promise[Unit])
}

class SimpleWorkflowActorSpec extends CromwellTestkitSpec {

  private def buildWorkflowActor(sampleWdl: SampleWdl, rawInputsOverride: String,
                                 workflowId: WorkflowId = WorkflowId.randomId(),
                                 matcher: Option[FailureMatcher] = None): WorkflowActorAndMetadataPromise = {
    val workflowSources = WorkflowSourceFiles(sampleWdl.wdlSource(), rawInputsOverride, "{}")
    val promise = Promise[Unit]()
    val watchActor = system.actorOf(MetadataWatchActor.props(matcher, promise), s"service-registry-$workflowId-${UUID.randomUUID()}")
    val workflowActor = TestFSMRef(new WorkflowActor(workflowId, StartNewWorkflow, workflowSources, ConfigFactory.load(),
      watchActor,
      system.actorOf(Props.empty, s"workflow-copy-log-router-$workflowId-${UUID.randomUUID()}"),
      system.actorOf(JobStoreActor.props(WriteCountingJobStoreDatabase.makeNew))),
      name = s"workflow-actor-$workflowId")
    WorkflowActorAndMetadataPromise(workflowActor, promise)
  }

  val TestExecutionTimeout = 10.seconds.dilated

  "A WorkflowActor" should {
    "start, run, succeed and die" in {
      startingCallsFilter("hello.hello") {
        val WorkflowActorAndMetadataPromise(workflowActor, _) = buildWorkflowActor(SampleWdl.HelloWorld, SampleWdl.HelloWorld.wdlJson)
        val probe = TestProbe()
        probe watch workflowActor
        within(TestExecutionTimeout) {
          waitForPattern("transition from ReadyToMaterializeState to MaterializationSuccessfulState") {
            waitForPattern("transition from FinalizingWorkflowState to WorkflowSucceededState") {
              workflowActor ! StartWorkflowCommand
            }
          }
        }
        probe.expectTerminated(workflowActor, 10.seconds.dilated)
      }
    }

    "fail to construct with missing inputs" in {
      val workflowId = WorkflowId.randomId()
      val failureMatcher = FailureMatcher("Required workflow input 'hello.hello.addressee' not specified.")
      val WorkflowActorAndMetadataPromise(workflowActor, promise) = buildWorkflowActor(SampleWdl.HelloWorld, "{}", workflowId, Option(failureMatcher))
      val probe = TestProbe()
      probe watch workflowActor
      within(TestExecutionTimeout) {
        waitForPattern("transition from MaterializingWorkflowDescriptorState to WorkflowFailedState") {
          workflowActor ! StartWorkflowCommand
          Await.result(promise.future, Duration.Inf)
        }
      }
      probe.expectTerminated(workflowActor, 10.seconds.dilated)
    }

    "fail to construct with inputs of the wrong type" in {
      val failureMatcher = FailureMatcher("Could not coerce value for 'hello.hello.addressee' into: WdlStringType")
      val WorkflowActorAndMetadataPromise(workflowActor, promise) = buildWorkflowActor(SampleWdl.HelloWorld, s""" { "$Addressee" : 3} """,
        WorkflowId.randomId(), Option(failureMatcher))

      val probe = TestProbe()
      probe watch workflowActor
      within(TestExecutionTimeout) {
        waitForPattern("transition from MaterializingWorkflowDescriptorState to WorkflowFailedState") {
          workflowActor ! StartWorkflowCommand
          Await.result(promise.future, Duration.Inf)
        }
      }
      probe.expectTerminated(workflowActor, 10.seconds.dilated)
    }

    "fail when a call fails" in {
      startingCallsFilter("goodbye.goodbye") {
        val WorkflowActorAndMetadataPromise(workflowActor, _) = buildWorkflowActor(SampleWdl.GoodbyeWorld, SampleWdl.GoodbyeWorld.wdlJson)
        val probe = TestProbe()
        probe watch workflowActor
        within(TestExecutionTimeout) {
          waitForPattern("transitioning from WorkflowExecutionInProgressState to WorkflowExecutionFailedState.") {
            waitForPattern("transitioning from FinalizationInProgressState to FinalizationSucceededState.") {
              workflowActor ! StartWorkflowCommand
            }
          }
        }
        probe.expectTerminated(workflowActor, 10.seconds.dilated)
      }
    }

    "gracefully handle malformed WDL" in {
      val WorkflowActorAndMetadataPromise(workflowActor, _) = buildWorkflowActor(SampleWdl.CoercionNotDefined, SampleWdl.CoercionNotDefined.wdlJson)
      val probe = TestProbe()
      probe watch workflowActor
      within(TestExecutionTimeout) {
        waitForPattern("transitioning from WorkflowExecutionInProgressState to WorkflowExecutionFailedState") {
          workflowActor ! StartWorkflowCommand
        }
      }
      probe.expectTerminated(workflowActor, 10.seconds.dilated)
    }
  }
}
