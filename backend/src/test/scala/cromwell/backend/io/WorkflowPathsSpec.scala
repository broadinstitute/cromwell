package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobBreadCrumb, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{JobKey, WorkflowId}
import cromwell.util.WomMocks
import org.mockito.Mockito._
import org.scalatest.{FlatSpec, Matchers}
import wom.graph.WomIdentifier

class WorkflowPathsSpec extends FlatSpec with Matchers with BackendSpec {

  def createBackendConfig(values: Map[String, String]): Config = {
    val backendConfig = mock[Config]
    values.foreach {
      case (key: String, value: String) => {
        when(backendConfig.hasPath(key)).thenReturn(true)
        when(backendConfig.getString(key)).thenReturn(value)
      }
    }
    backendConfig
  }

  def rootConfig(root: Option[String], dockerRoot: Option[String]) = {
    val values: Map[String,String] = root.map("root" -> _).toMap ++ dockerRoot.map("dockerRoot" -> _).toMap
    createBackendConfig(values)
  }

  def testWorkflowPaths(root: Option[String], dockerRoot: Option[String]): Unit = {

  }

  "WorkflowPaths" should "provide correct paths for a workflow" in {
    val backendConfig = createBackendConfig(Map("root" -> "local-cromwell-executions"))
    val wd = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = new WorkflowPathsWithDocker(wd, backendConfig)
    val id = wd.id
    workflowPaths.workflowRoot.pathAsString shouldBe
      DefaultPathBuilder.get(s"local-cromwell-executions/wf_hello/$id").toAbsolutePath.pathAsString
    workflowPaths.dockerWorkflowRoot.pathAsString shouldBe
      s"/cromwell-executions/wf_hello/$id"
  }

  "WorkflowPaths" should "provide correct paths for a sub workflow" in {
    val backendConfig = createBackendConfig(Map("root" -> "local-cromwell-executions"))
    when(backendConfig.hasPath("root")).thenReturn(true)
    when(backendConfig.getString("root")).thenReturn("local-cromwell-executions") // This is the folder defined in the config as the execution root dir
    val rootWd = mock[BackendWorkflowDescriptor]
    val rootWorkflow = WomMocks.mockWorkflowDefinition("rootWorkflow")
    val rootWorkflowId = WorkflowId.randomId()
    rootWd.callable returns rootWorkflow
    rootWd.id returns rootWorkflowId

    val subWd = mock[BackendWorkflowDescriptor]
    val subWorkflow = WomMocks.mockWorkflowDefinition("subWorkflow")
    val subWorkflowId = WorkflowId.randomId()
    subWd.callable returns subWorkflow
    subWd.id returns subWorkflowId
    
    val call1 = WomMocks.mockTaskCall(WomIdentifier("call1"))
    
    val jobKey = new JobKey {
      override def node = call1
      override def tag: String = "tag1"
      override def index: Option[Int] = Option(1)
      override def attempt: Int = 2
    }

    subWd.breadCrumbs returns List(BackendJobBreadCrumb(rootWorkflow, rootWorkflowId, jobKey))
    subWd.id returns subWorkflowId
    
    val workflowPaths = new WorkflowPathsWithDocker(subWd, backendConfig)
    workflowPaths.workflowRoot.pathAsString shouldBe
      DefaultPathBuilder.get(
        s"local-cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
      ).toAbsolutePath.pathAsString
    workflowPaths.dockerWorkflowRoot.pathAsString shouldBe s"/cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
  }
}
