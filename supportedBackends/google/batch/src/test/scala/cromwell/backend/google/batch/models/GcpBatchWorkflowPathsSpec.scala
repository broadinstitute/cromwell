package cromwell.backend.google.batch.models

import com.google.cloud.NoCredentials
import common.collections.EnhancedCollections._
import cromwell.backend.google.batch.actors.GcpBatchInitializationActor
import cromwell.backend.{BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsString}

class GcpBatchWorkflowPathsSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {
  import BackendSpec._
  import GcpBatchTestConfig._
  import cromwell.filesystems.gcs.MockGcsPathBuilder._

  behavior of "GcpBatchWorkflowPaths"

  var workflowDescriptor: BackendWorkflowDescriptor = _
  var workflowPaths: GcpBatchWorkflowPaths = _

  override def beforeAll(): Unit = {
    workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    workflowPaths = GcpBatchWorkflowPaths(
      workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      gcpBatchConfiguration,
      pathBuilders(),
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )
  }

  it should "map the correct paths" in {
    workflowPaths.executionRoot.pathAsString should be("gs://my-cromwell-workflows-bucket/")
    workflowPaths.workflowRoot.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/")
  }
}
