package cromwell.backend.io

import better.files._
import com.typesafe.config.Config
import cromwell.backend.{BackendJobBreadCrumb, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.{JobKey, WorkflowId}
import org.mockito.Mockito._
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.{Call, Workflow}

class WorkflowPathsSpec extends FlatSpec with Matchers with BackendSpec {

  val backendConfig = mock[Config]

  "WorkflowPaths" should "provide correct paths for a workflow" in {
    when(backendConfig.hasPath(any[String])).thenReturn(true)
    when(backendConfig.getString(any[String])).thenReturn("local-cromwell-executions") // This is the folder defined in the config as the execution root dir
    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = new WorkflowPathsWithDocker(wd, backendConfig)
    val id = wd.id
    workflowPaths.workflowRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id").pathAsString
    workflowPaths.dockerWorkflowRoot.toString shouldBe
      s"/root/wf_hello/$id"
  }

  "WorkflowPaths" should "provide correct paths for a sub workflow" in {
    when(backendConfig.hasPath(any[String])).thenReturn(true)
    when(backendConfig.getString(any[String])).thenReturn("local-cromwell-executions") // This is the folder defined in the config as the execution root dir
    
    val rootWd = mock[BackendWorkflowDescriptor]
    val rootWorkflow = mock[Workflow]
    val rootWorkflowId = WorkflowId.randomId()
    rootWorkflow.unqualifiedName returns "rootWorkflow"
    rootWd.workflow returns rootWorkflow
    rootWd.id returns rootWorkflowId

    val subWd = mock[BackendWorkflowDescriptor]
    val subWorkflow = mock[Workflow]
    val subWorkflowId = WorkflowId.randomId()
    subWorkflow.unqualifiedName returns "subWorkflow"
    subWd.workflow returns subWorkflow
    subWd.id returns subWorkflowId
    
    val call1 = mock[Call]
    call1.unqualifiedName returns "call1"
    val call2 = mock[Call]
    call2.unqualifiedName returns "call2"
    
    val jobKey = new JobKey {
      override def scope = call1
      override def tag: String = "tag1"
      override def index: Option[Int] = Option(1)
      override def attempt: Int = 2
    }

    subWd.breadCrumbs returns List(BackendJobBreadCrumb(rootWorkflow, rootWorkflowId, jobKey))
    subWd.id returns subWorkflowId
    
    val workflowPaths = new WorkflowPathsWithDocker(subWd, backendConfig)
    workflowPaths.workflowRoot.toString shouldBe File(s"local-cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId").pathAsString
    workflowPaths.dockerWorkflowRoot.toString shouldBe s"/root/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
  }
}
