package cromwell.engine.workflow.lifecycle

import akka.testkit.EventFilter
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends, WorkflowDescriptorBuilder}
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor.StartExecutingWorkflowCommand
import cromwell.util.SampleWdl
import org.scalatest.BeforeAndAfter

class WorkflowExecutionActorSpec extends CromwellTestkitSpec with BeforeAndAfter with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  val stubbedConfig = ConfigFactory.load().getConfig("backend.providers.Mock")

  val MockBackendConfigEntry = BackendConfigurationEntry(
    name = "Mock",
    lifecycleActorFactoryClass = "cromwell.engine.backend.mock.RetryableBackendLifecycleActorFactory",
    stubbedConfig.getConfig("config")
  )

  CromwellBackends.initBackends(List(MockBackendConfigEntry), MockBackendConfigEntry, system)

  val runtimeSection =
    """
      |runtime {
      | backend: "Mock"
      |}
    """.stripMargin

  val workflowId = WorkflowId.randomId()
  val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.HelloWorld.asWorkflowSources(runtime = runtimeSection))

  "WorkflowExecutionActor" should {
    "retry a job 2 times and succeed in the third attempt" in {
      val workflowExecutionActor = system.actorOf(WorkflowExecutionActor.props(workflowId, engineWorkflowDescriptor), "WorkflowExecutionActor")
      EventFilter.info(pattern = ".*Final Outputs", occurrences = 1).intercept {
        EventFilter.info(pattern = "Starting calls: hello.hello", occurrences = 3).intercept {
          workflowExecutionActor ! StartExecutingWorkflowCommand
        }
      }
    }
  }

  override def afterAll(): Unit = {
    system.shutdown()
  }
}
