package cromwell.engine.workflow.lifecycle

import akka.actor.Actor
import akka.testkit.{TestActorRef, EventFilter}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends}
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.StartExecutingWorkflowCommand
import cromwell.util.SampleWdl
import org.scalatest.BeforeAndAfter

class WorkflowExecutionActorSpec extends CromwellTestkitSpec with BeforeAndAfter with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  def mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  val stubbedConfig = ConfigFactory.load().getConfig("backend.providers.Mock").getConfig("config")

  val runtimeSection =
    """
      |runtime {
      | backend: "Mock"
      |}
    """.stripMargin

  "WorkflowExecutionActor" should {
    "retry a job 2 times and succeed in the third attempt" in {
      val serviceRegistry = mockServiceRegistryActor
      val MockBackendConfigEntry = BackendConfigurationEntry(
        name = "Mock",
        lifecycleActorFactoryClass = "cromwell.engine.backend.mock.RetryableBackendLifecycleActorFactory",
        stubbedConfig
      )
      CromwellBackends.initBackends(List(MockBackendConfigEntry), MockBackendConfigEntry, system)

      val workflowId = WorkflowId.randomId()
      val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.HelloWorld.asWorkflowSources(runtime = runtimeSection))
      val workflowExecutionActor = system.actorOf(WorkflowExecutionActor.props(workflowId, engineWorkflowDescriptor, serviceRegistry), "WorkflowExecutionActor")

      EventFilter.info(pattern = ".*Final Outputs", occurrences = 1).intercept {
        EventFilter.info(pattern = "Starting calls: hello.hello", occurrences = 3).intercept {
          workflowExecutionActor ! StartExecutingWorkflowCommand
        }
      }
      system.stop(serviceRegistry)
    }

    "execute a workflow with scatters" in {
      val serviceRegistry = mockServiceRegistryActor
      val MockBackendConfigEntry = BackendConfigurationEntry(
        name = "Mock",
        lifecycleActorFactoryClass = "cromwell.engine.backend.mock.DefaultBackendLifecycleActorFactory",
        stubbedConfig
      )
      CromwellBackends.initBackends(List(MockBackendConfigEntry), MockBackendConfigEntry, system)

      val workflowId = WorkflowId.randomId()
      val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.SimpleScatterWdl.asWorkflowSources(runtime = runtimeSection))
      val workflowExecutionActor = system.actorOf(WorkflowExecutionActor.props(workflowId, engineWorkflowDescriptor, serviceRegistry), "WorkflowExecutionActor")

      val scatterLog = "Starting calls: scatter0.inside_scatter:0:1, scatter0.inside_scatter:1:1, scatter0.inside_scatter:2:1, scatter0.inside_scatter:3:1, scatter0.inside_scatter:4:1"

      EventFilter.info(pattern = ".*Final Outputs", occurrences = 1).intercept {
        EventFilter.info(pattern = scatterLog, occurrences = 1).intercept {
          EventFilter.info(pattern = "Starting calls: scatter0.outside_scatter:NA:1", occurrences = 1).intercept {
            workflowExecutionActor ! StartExecutingWorkflowCommand
          }
        }
      }
      system.stop(serviceRegistry)
    }
  }

  override def afterAll(): Unit = {
    system.shutdown()
  }
}
