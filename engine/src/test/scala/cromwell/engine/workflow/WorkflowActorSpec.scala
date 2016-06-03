package cromwell.engine.workflow

import akka.actor.Actor
import akka.testkit.{TestActorRef, TestFSMRef}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.{CromwellBackends, WorkflowDescriptorBuilder}
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{WorkflowDescriptorMaterializationResult, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.WorkflowInitializationSucceededResponse
import cromwell.util.SampleWdl.{HelloWorld, ThreeStep}
import org.scalatest.BeforeAndAfter
import spray.json.DefaultJsonProtocol._
import spray.json._


class WorkflowActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder with BeforeAndAfter {
  override implicit val actorSystem = system

  val mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  private def createWorkflowActor(workflowId: WorkflowId = WorkflowId.randomId(),
                                  startMode: StartMode = StartNewWorkflow,
                                  wdlSource: WorkflowSourceFiles,
                                  messageStopper: PartialFunction[Any, Unit]) = {
    TestFSMRef(new WorkflowActor(workflowId, startMode, wdlSource, ConfigFactory.load, mockServiceRegistryActor) {
      override def receive = messageStopper orElse super.receive
    })
  }

  before {
    val local = CromwellTestkitSpec.DefaultLocalBackendConfigEntry
    CromwellBackends.initBackends(List(local), local, system)
  }

  "ShadowWorkflowActor" should {
    "move from WorkflowUnstartedState to MaterializingWorkflowDescriptorState" in {
      val wdlSources = HelloWorld.asWorkflowSources()
      val messageStopper: PartialFunction[Any, Unit] = {
        case _: WorkflowDescriptorMaterializationResult =>
      }
      val testActor = createWorkflowActor(wdlSource = wdlSources, messageStopper = messageStopper)
      testActor.setState(stateName = WorkflowUnstartedState, stateData = WorkflowActorData.empty)
      testActor ! StartWorkflowCommand
      testActor.stateName should be(WorkflowActor.MaterializingWorkflowDescriptorState)
      testActor.stop()
    }

    "transition to InitializingWorkflowState" in {
      import scala.concurrent.duration._

      // WorkflowOptions (local) should trump runtime attributes (SGE)
      val wdlSources = ThreeStep.asWorkflowSources()
      val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = wdlSources)
      val messageStopper: PartialFunction[Any, Unit] = {
        case WorkflowInitializationSucceededResponse =>
      }
      val testActor = createWorkflowActor(wdlSource = wdlSources, messageStopper = messageStopper)
      testActor.setState(stateName = MaterializingWorkflowDescriptorState, stateData = WorkflowActorData.empty)
      within(5.seconds) {
        testActor ! MaterializeWorkflowDescriptorSuccessResponse(descriptor)
        testActor.stateName should be(WorkflowActor.InitializingWorkflowState)
        testActor.stateData.workflowDescriptor should be(Some(descriptor))
        testActor.stateData.currentLifecycleStateActor.isDefined should be(true)
      }
      testActor.stop()
    }
  }
}
