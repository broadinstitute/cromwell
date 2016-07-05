package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.BackendSpec
import org.mockito.Matchers._
import org.mockito.Mockito._
import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

class WorkflowPathsSpec extends FlatSpec with Matchers with BackendSpec with MockitoSugar {

  val backendConfig = mock[Config]

  "WorkflowPaths" should "provide correct paths for a workflow" in {
    when(backendConfig.getString(any[String])).thenReturn("local-cromwell-executions") // This is the folder defined in the config as the execution root dir
    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = new WorkflowPaths(wd, backendConfig, None)
    val id = wd.id
    workflowPaths.workflowRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id"
    workflowPaths.dockerWorkflowRoot.toString shouldBe
      s"/root/hello/$id"
  }
}
