package cromwell.filesystems.blob

import com.typesafe.config.ConfigFactory
import common.exception.AggregatedMessageException
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class BlobFileSystemConfigSpec extends AnyFlatSpec with Matchers {

  private val workspaceManagerURL = WorkspaceManagerURL("https://wsm.example.com")
  private val b2cToken = "b0gus-t0ken"

  it should "parse configs for a minimal functioning factory with native blob access" in {
    val config = BlobFileSystemConfig(
      ConfigFactory.parseString(s"""
      """.stripMargin)
    )
    config.expiryBufferMinutes should equal(BlobFileSystemConfig.defaultExpiryBufferMinutes)
  }

  it should "parse configs for a functioning factory with WSM-mediated blob access" in {
    val config = BlobFileSystemConfig(
      ConfigFactory.parseString(s"""
                                   |expiry-buffer-minutes = "20"
                                   |workspace-manager {
                                   |  url = "$workspaceManagerURL"
                                   |  b2cToken = "$b2cToken"
                                   |}
                                   |
        """.stripMargin)
    )
    config.expiryBufferMinutes should equal(20L)
    config.workspaceManagerConfig.isDefined shouldBe true
    config.workspaceManagerConfig.get.url shouldBe workspaceManagerURL
    config.workspaceManagerConfig.get.overrideWsmAuthToken.contains(b2cToken) shouldBe true
  }

  it should "fail when partial WSM config is supplied" in {
    val rawConfig =
      ConfigFactory.parseString(s"""
                                   |expiry-buffer-minutes = "10"
                                   |workspace-manager {
                                   |  b2cToken = "$b2cToken"
                                   |}
                                   |
        """.stripMargin)

    val error = intercept[AggregatedMessageException](BlobFileSystemConfig(rawConfig))
    error.getMessage should include("No configuration setting found for key 'url'")
  }
}
