package cromwell.engine.workflow

import akka.actor.Actor
import akka.testkit.{TestActorRef, TestFSMRef}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.{CromwellBackends, WorkflowDescriptorBuilder}
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.MaterializeWorkflowDescriptorSuccessResponse
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
                                  wdlSource: WorkflowSourceFiles) = {
    TestFSMRef(new WorkflowActor(workflowId, startMode, wdlSource, ConfigFactory.load, mockServiceRegistryActor))
  }

  before {
    val local = CromwellTestkitSpec.DefaultLocalBackendConfigEntry
    CromwellBackends.initBackends(List(local), local, system)
  }

  "ShadowWorkflowActor" should {
    "move from WorkflowUnstartedState to MaterializingWorkflowDescriptorState" in {
      val wdlSources = WorkflowSourceFiles(HelloWorld.wdlSource("runtime { }"),
        HelloWorld.rawInputs.toJson.toString(), "{ }")
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = WorkflowUnstartedState, stateData = WorkflowActorData.empty)
      testActor ! StartWorkflowCommand
      testActor.stateName should be (WorkflowActor.MaterializingWorkflowDescriptorState)
      testActor.stop()
    }

    "transition to InitializingWorkflowState with correct call assignments given workflow options" in {
      import scala.concurrent.duration._

      val wfOptions =
        """{
          | "backend": "local"
          |}
        """.stripMargin
      val runtimeAttributes =
        """runtime {
          | backend: "SGE"
          |}""".stripMargin
      // WorkflowOptions (local) should trump runtime attributes (SGE)
      val wdlSources = ThreeStep.asWorkflowSources(runtime = runtimeAttributes, workflowOptions = wfOptions)
      val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = wdlSources)
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = MaterializingWorkflowDescriptorState, stateData = WorkflowActorData.empty)
      within(5.seconds) {
        testActor ! MaterializeWorkflowDescriptorSuccessResponse(descriptor)
        testActor.stateName should be(WorkflowActor.InitializingWorkflowState)
        testActor.stateData.workflowDescriptor should be(Some(descriptor))
        testActor.stateData.currentLifecycleStateActor.isDefined should be(true)
      }
      testActor.stop()
    }

    "transition to InitializingWorkflowState with correct call assignments given runtime-attributes" in {
      import scala.concurrent.duration._
      val runtimeAttributes =
        """runtime {
          | backend: "local"
          |}""".stripMargin
      val wdlSources = ThreeStep.asWorkflowSources(runtime = runtimeAttributes)
      val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = wdlSources)
      val testActor = createWorkflowActor(wdlSource = wdlSources)
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
