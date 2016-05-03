package cromwell.engine.workflow.lifecycle

import akka.testkit.EventFilter
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends, WorkflowDescriptorBuilder}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.StartExecutingWorkflowCommand
import cromwell.util.SampleWdl
import org.scalatest.BeforeAndAfter

class WorkflowExecutionActorSpec extends CromwellTestkitSpec with BeforeAndAfter with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  val stubbedConfig = ConfigFactory.load().getConfig("backend.providers.Mock").getConfig("config")

  val runtimeSection =
    """
      |runtime {
      | backend: "Mock"
      |}
    """.stripMargin

  "WorkflowExecutionActor" should {
    "retry a job 2 times and succeed in the third attempt" in {
      val MockBackendConfigEntry = BackendConfigurationEntry(
        name = "Mock",
        lifecycleActorFactoryClass = "cromwell.engine.backend.mock.RetryableBackendLifecycleActorFactory",
        stubbedConfig
      )
      CromwellBackends.initBackends(List(MockBackendConfigEntry), MockBackendConfigEntry, system)

      val workflowId = WorkflowId.randomId()
      val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.HelloWorld.asWorkflowSources(runtime = runtimeSection))
      val workflowExecutionActor = system.actorOf(WorkflowExecutionActor.props(workflowId, engineWorkflowDescriptor), "WorkflowExecutionActor")

      EventFilter.info(pattern = ".*Final Outputs", occurrences = 1).intercept {
        EventFilter.info(pattern = "Starting calls: hello.hello", occurrences = 3).intercept {
          workflowExecutionActor ! StartExecutingWorkflowCommand
        }
      }
    }

    "execute a workflow with scatters" in {
      val MockBackendConfigEntry = BackendConfigurationEntry(
        name = "Mock",
        lifecycleActorFactoryClass = "cromwell.engine.backend.mock.DefaultBackendLifecycleActorFactory",
        stubbedConfig
      )
      CromwellBackends.initBackends(List(MockBackendConfigEntry), MockBackendConfigEntry, system)

      val workflowId = WorkflowId.randomId()
      val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.SimpleScatterWdl.asWorkflowSources(runtime = runtimeSection))
      val workflowExecutionActor = system.actorOf(WorkflowExecutionActor.props(workflowId, engineWorkflowDescriptor), "WorkflowExecutionActor")

      val scatterLog = "Starting calls: scatter0.inside_scatter:0:1, scatter0.inside_scatter:1:1, scatter0.inside_scatter:2:1, scatter0.inside_scatter:3:1, scatter0.inside_scatter:4:1"

      EventFilter.info(pattern = ".*Final Outputs", occurrences = 1).intercept {
        EventFilter.info(pattern = scatterLog, occurrences = 1).intercept {
          EventFilter.info(pattern = "Starting calls: scatter0.outside_scatter:NA:1", occurrences = 1).intercept {
            workflowExecutionActor ! StartExecutingWorkflowCommand
          }
        }
      }
    }
  }

  override def afterAll(): Unit = {
    system.shutdown()
  }
}
