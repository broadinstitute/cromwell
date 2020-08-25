package cromwell.backend.google.pipelines.common

import com.google.cloud.NoCredentials
import common.collections.EnhancedCollections._
import cromwell.backend.{BackendSpec, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.auth.GoogleAuthModeSpec
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsString}

class PipelinesApiWorkflowPathsSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Mockito {
  import BackendSpec._
  import PipelinesApiTestConfig._
  import cromwell.filesystems.gcs.MockGcsPathBuilder._

  behavior of "PipelinesApiWorkflowPaths"

  var workflowDescriptor: BackendWorkflowDescriptor = _
  var workflowPaths: PipelinesApiWorkflowPaths = _

  override def beforeAll: Unit = {
    workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    workflowPaths = PipelinesApiWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), papiConfiguration, pathBuilders, PipelinesApiInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper)
  }

  it should "map the correct paths" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    workflowPaths.executionRoot.pathAsString should be("gs://my-cromwell-workflows-bucket/")
    workflowPaths.workflowRoot.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/")
    workflowPaths.gcsAuthFilePath.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/${workflowDescriptor.id}_auth.json")
  }

  it should "calculate the call cache path prefix from the workflow execution root correctly" in {
    val WorkspaceBucket = "gs://workspace-id"
    val ExecutionRoot = WorkspaceBucket + "/submission-id"
    PipelinesApiWorkflowPaths.callCachePathPrefixFromExecutionRoot(ExecutionRoot) shouldBe WorkspaceBucket
  }
}
