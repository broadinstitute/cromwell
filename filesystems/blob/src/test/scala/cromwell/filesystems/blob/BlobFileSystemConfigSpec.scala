package cromwell.filesystems.blob

import com.typesafe.config.ConfigFactory
import common.exception.AggregatedMessageException
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.util.UUID

class BlobFileSystemConfigSpec extends AnyFlatSpec with Matchers {

  private val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
  private val container = BlobContainerName("storageContainer")
  private val workspaceId = WorkspaceId(UUID.fromString("B0BAFE77-0000-0000-0000-000000000000"))
  private val containerResourceId = ContainerResourceId(UUID.fromString("F00B4911-0000-0000-0000-000000000000"))
  private val workspaceManagerURL = WorkspaceManagerURL("https://wsm.example.com")
  private val b2cToken = "b0gus-t0ken"

  it should "parse configs for a minimal functioning factory with native blob access" in {
    val config = BlobFileSystemConfig(
      ConfigFactory.parseString(
      s"""
         |container = "$container"
         |endpoint = "$endpoint"
      """.stripMargin)
    )
    config.blobContainerName should equal(container)
    config.endpointURL should equal(endpoint)
    config.expiryBufferMinutes should equal(BlobFileSystemConfig.defaultExpiryBufferMinutes)
  }

  it should "parse configs for a functioning factory with WSM-mediated blob access" in {
    val config = BlobFileSystemConfig(
      ConfigFactory.parseString(
        s"""
           |container = "$container"
           |endpoint = "$endpoint"
           |expiry-buffer-minutes = "20"
           |workspace-manager {
           |  url = "$workspaceManagerURL"
           |  workspace-id = "$workspaceId"
           |  container-resource-id = "$containerResourceId"
           |  b2cToken = "$b2cToken"
           |}
           |
        """.stripMargin)
    )
    config.blobContainerName should equal(container)
    config.endpointURL should equal(endpoint)
    config.expiryBufferMinutes should equal(20L)
    config.workspaceManagerConfig.isDefined shouldBe true
    config.workspaceManagerConfig.get.url shouldBe workspaceManagerURL
    config.workspaceManagerConfig.get.workspaceId shouldBe workspaceId
    config.workspaceManagerConfig.get.containerResourceId shouldBe containerResourceId
    config.workspaceManagerConfig.get.overrideWsmAuthToken.contains(b2cToken) shouldBe true
  }

  it should "fail when partial WSM config is supplied" in {
    val rawConfig =
      ConfigFactory.parseString(
        s"""
         |container = "$container"
         |endpoint = "$endpoint"
         |expiry-buffer-minutes = "10"
         |workspace-manager {
         |  url = "$workspaceManagerURL"
         |  container-resource-id = "$containerResourceId"
         |}
         |
        """.stripMargin)

    val error = intercept[AggregatedMessageException](BlobFileSystemConfig(rawConfig))
    error.getMessage should include("No configuration setting found for key 'workspace-id'")
  }
}
