package cromwell.backend.impl.jes

import cromwell.backend.BackendSpec
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import cromwell.backend.impl.jes.MockObjects._

class JesWorkflowPathsSpec extends FlatSpec with Matchers with Mockito {
  import BackendSpec._
  import JesTestConfig._

  behavior of "JesWorkflowPaths"

  it should "map the correct paths" in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val workflowPaths = JesWorkflowPaths(workflowDescriptor, jesConfiguration)
    workflowPaths.rootPath.toUri.toString should be("gs://my-cromwell-workflows-bucket/")
    workflowPaths.workflowRootPath.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/")
    workflowPaths.gcsAuthFilePath.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/${workflowDescriptor.id}_auth.json")
  }
}
