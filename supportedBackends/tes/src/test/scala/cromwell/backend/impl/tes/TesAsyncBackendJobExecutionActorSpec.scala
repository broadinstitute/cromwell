package cromwell.backend.impl.tes

import com.fasterxml.jackson.databind.ext.NioPathSerializer
import cromwell.core.path.{NioPath, Path}
import cromwell.filesystems.blob.{BlobPath, BlobPathBuilderFactory, EndpointURL}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.tools.nsc.io.Path


class TesAsyncBackendJobExecutionActorSpec extends AnyFlatSpec with Matchers {
  behavior of "TesAsyncBackendJobExecutionActor"

  val fullyQualifiedName = "this.name.is.more.than.qualified"
  val workflowName = "mockWorkflow"
  val someBlobUrl = "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-d8143fd8-aa07-446d-9ba0-af72203f1794/nyxp6c/tes-internal/configuration/supported-vm-sizes?sv=sasToken"
  val someNotBlobUrl = "https://www.google.com/path/to/exile"
  var index = 0

  val blobInput_0 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someBlobUrl),
    path = someBlobUrl,
    `type` = Option("FILE"),
    content = None
  )
  index = index+1

  val blobInput_1 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someBlobUrl),
    path = someBlobUrl,
    `type` = Option("FILE"),
    content = None
  )
  index = index+1

  val notBlobInput_1 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someNotBlobUrl + index),
    path = someNotBlobUrl + index,
    `type` = Option("FILE"),
    content = None
  )
  index = index+1

  val notBlobInput_2 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someNotBlobUrl + index),
    path = someNotBlobUrl + index,
    `type` = Option("FILE"),
    content = None
  )

  def pathGetter(pathString: String): Path = {
    /*
    val factory: BlobPathBuilderFactory = BlobPathBuilderFactory()
    val nioPath: java.nio.file.Path = NioPathSerializer()
    val endpoint: EndpointURL = EndpointURL("www.some.blob.endpoint/url")
    if(pathString.startsWith(someBlobUrl)) BlobPath("someBlobPath",endpoint,)
     */
  }

  it should "not generate sas params when no blob paths are provided" in {
    val inputs = List[Input]
    //TesAsyncBackendJobExecutionActor.getLocalizedSasTokenParams(inputs,)
  }

  it should "generate proper script preamble" in {
    val mockEndpoint = "www.workspacemanager.com"
    val mockWorkspaceId = "1111-2222-3333-4444"
    val mockContainerId = "5678-who-do-we-appreciate"
    val mockSasParams: LocalizedSasTokenParams = LocalizedSasTokenParams(mockEndpoint, mockWorkspaceId, mockContainerId)
    TesAsyncBackendJobExecutionActor.generateLocalizedSasScriptPreammble(mockSasParams) shouldBe
      s"""
         |WSM_ENDPOINT="${mockEndpoint}"
         |WORKSPACE_ID="${mockWorkspaceId}"
         |CONTAINER_RESOURCE_ID="${mockContainerId}"
         |""".stripMargin
  }
}