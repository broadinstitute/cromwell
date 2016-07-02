package cromwell.backend.impl.jes

import cromwell.backend.BackendSpec
import cromwell.core.IntegrationTest
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}

class JesWorkflowPathsSpec extends FlatSpec with Matchers {
  import BackendSpec._
  import JesTestConfig._

  behavior of "JesWorkflowPaths"

  it should "map the correct paths" taggedAs IntegrationTest in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val workflowPaths = JesWorkflowPaths(workflowDescriptor, jesConfiguration)
    workflowPaths.rootPath.toString should be("gs://my-cromwell-workflows-bucket")
    workflowPaths.workflowRootPath.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}")
    workflowPaths.gcsAuthFilePath.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/${workflowDescriptor.id}_auth.json")
  }
}
