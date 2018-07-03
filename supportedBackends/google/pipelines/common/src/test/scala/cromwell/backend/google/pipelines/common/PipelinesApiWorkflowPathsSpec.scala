package cromwell.backend.google.pipelines.common

import com.google.cloud.NoCredentials
import cromwell.backend.BackendSpec
import cromwell.cloudsupport.gcp.auth.GoogleAuthModeSpec
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsString}

class PipelinesApiWorkflowPathsSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {
  import BackendSpec._
  import PipelinesApiTestConfig._
  import cromwell.filesystems.gcs.MockGcsPathBuilder._

  behavior of "JesWorkflowPaths"

  it should "map the correct paths" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.mapValues(JsString.apply)).compactPrint)
    )
    val workflowPaths = PipelinesApiWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration, pathBuilders)
    workflowPaths.executionRoot.pathAsString should be("gs://my-cromwell-workflows-bucket/")
    workflowPaths.workflowRoot.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/")
    workflowPaths.gcsAuthFilePath.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/${workflowDescriptor.id}_auth.json")
  }
}
