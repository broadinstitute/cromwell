package cromwell.backend.impl.local

import org.scalatest.{FlatSpec, Matchers}

class WorkflowPathsSpec extends FlatSpec with Matchers with BackendTestkitSpec {

  "WorkflowPaths" should "provide correct paths for a workflow" in {

    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = new WorkflowPaths(wd, defaultConfig)
    val id = wd.id
    workflowPaths.workflowRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id"
    workflowPaths.dockerWorkflowRoot.toString shouldBe
      s"/root/hello/$id"
  }
}
