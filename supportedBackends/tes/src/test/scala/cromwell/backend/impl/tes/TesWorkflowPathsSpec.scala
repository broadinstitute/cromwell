package cromwell.backend.impl.tes

import better.files._
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.{BackendJobBreadCrumb, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.{JobKey, WorkflowId}
import cromwell.util.WomMocks
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.graph.WomIdentifier

class TesWorkflowPathsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with BackendSpec {

  "WorkflowPaths" should "provide correct paths for a workflow" in {
    val wd = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = TesWorkflowPaths(wd, TesTestConfig.backendConfig)
    val id = wd.id
    workflowPaths.workflowRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id").pathAsString
    workflowPaths.dockerWorkflowRoot.toString shouldBe
      s"/cromwell-executions/wf_hello/$id"
  }

  "WorkflowPaths" should "provide correct paths for a sub workflow" in {
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
    
    val workflowPaths = TesWorkflowPaths(subWd, TesTestConfig.backendConfig)
    workflowPaths.workflowRoot.toString shouldBe File(s"local-cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId").pathAsString
    workflowPaths.dockerWorkflowRoot.toString shouldBe s"/cromwell-executions/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
  }
}
