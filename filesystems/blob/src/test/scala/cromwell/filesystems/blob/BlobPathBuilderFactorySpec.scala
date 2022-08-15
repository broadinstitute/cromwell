package cromwell.filesystems.blob

import com.typesafe.config.ConfigFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class BlobPathBuilderFactorySpec extends AnyFlatSpec with Matchers {

  it should "parse configs for a functioning factory" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val store = "inputs"
    val sasToken = "{SAS TOKEN HERE}"
    val workspaceId = "mockWorkspaceId"
    val workspaceManagerURL = "https://test.ws.org"
    val instanceConfig = ConfigFactory.parseString(
      s"""
      |sas-token = "$sasToken"
      |store = "$store"
      |endpoint = "$endpoint"
      |workspace-id = "$workspaceId"
      """.stripMargin)
    val singletonConfig = ConfigFactory.parseString(s"""workspace-manager-url = "$workspaceManagerURL" """)
    val globalConfig = ConfigFactory.parseString("""""")
    val factory = BlobPathBuilderFactory(globalConfig, instanceConfig, new BlobFileSystemConfig(singletonConfig))
    factory.container should equal(store)
    factory.endpoint should equal(endpoint)
    factory.sasToken should equal(sasToken)
    factory.workspaceId should equal(Some(workspaceId))
    factory.workspaceManagerURL should equal(Some(workspaceManagerURL))
  }
}
