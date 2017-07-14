package cromwell.backend.impl.tes

import better.files._
import cromwell.backend.{BackendJobBreadCrumb, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.{JobKey, WorkflowId}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.{WdlCall, WdlWorkflow}

class TesWorkflowPathsSpec extends FlatSpec with Matchers with BackendSpec {

  "WorkflowPaths" should "provide correct paths for a workflow" in {
    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = TesWorkflowPaths(wd, TesTestConfig.backendConfig)
    val id = wd.id
    workflowPaths.workflowRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id").pathAsString
    workflowPaths.dockerWorkflowRoot.toString shouldBe
      s"/cromwell-executions/wf_hello/$id"
  }

  "WorkflowPaths" should "provide correct paths for a sub workflow" in {
    val rootWd = mock[BackendWorkflowDescriptor]
    val rootWorkflow = mock[WdlWorkflow]
    val rootWorkflowId = WorkflowId.randomId()
    rootWorkflow.unqualifiedName returns "rootWorkflow"
    rootWd.workflow returns rootWorkflow
    rootWd.id returns rootWorkflowId

    val subWd = mock[BackendWorkflowDescriptor]
    val subWorkflow = mock[WdlWorkflow]
    val subWorkflowId = WorkflowId.randomId()
    subWorkflow.unqualifiedName returns "subWorkflow"
    subWd.workflow returns subWorkflow
    subWd.id returns subWorkflowId
    
    val call1 = mock[WdlCall]
    call1.unqualifiedName returns "call1"
    val call2 = mock[WdlCall]
    call2.unqualifiedName returns "call2"
    
    val jobKey = new JobKey {
      override def scope = call1
      override def tag: String = "tag1"
      override def index: Option[Int] = Option(1)
      override def attempt: Int = 2
    }

    subWd.breadCrumbs returns List(BackendJobBreadCrumb(rootWorkflow, rootWorkflowId, jobKey))
    subWd.id returns subWorkflowId
    
    val workflowPaths = new TesWorkflowPaths(subWd, TesTestConfig.backendConfig)
    workflowPaths.workflowRoot.toString shouldBe File(s"local-cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId").pathAsString
    workflowPaths.dockerWorkflowRoot.toString shouldBe s"/cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
  }
}
