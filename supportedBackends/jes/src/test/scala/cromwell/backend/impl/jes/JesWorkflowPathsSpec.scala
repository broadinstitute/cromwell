package cromwell.backend.impl.jes

import cromwell.backend.BackendSpec
import cromwell.core.TestKitSuite
import cromwell.core.path.PathImplicits._
import cromwell.filesystems.gcs.auth.GoogleAuthModeSpec
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

class JesWorkflowPathsSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {
  import BackendSpec._
  import JesTestConfig._

  behavior of "JesWorkflowPaths"

  it should "map the correct paths" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val workflowPaths = JesWorkflowPaths(workflowDescriptor, jesConfiguration)(system)
    workflowPaths.executionRoot.toRealString should be("gs://my-cromwell-workflows-bucket/")
    workflowPaths.workflowRoot.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/")
    workflowPaths.gcsAuthFilePath.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/${workflowDescriptor.id}_auth.json")
  }
}
