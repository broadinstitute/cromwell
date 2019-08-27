package cromwell.filesystems.oss.nio

import java.net.URI
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.scalatest.mockito.MockitoSugar


object TTLOssStorageConfigurationSpec {

  val BcsBackendConfigString =
    s"""
     |  auth {
     |      endpoint = "oss-cn-shanghai.aliyuncs.com"
     |      access-id = "test-access-id"
     |      access-key = "test-access-key"
     |      security-token = "test-security-token"
     |  }
     |  caching {
     |      duplication-strategy = "reference"
     |  }
      """.stripMargin

  val BcsBackendConfig = ConfigFactory.parseString(BcsBackendConfigString)
}

class TTLOssStorageConfigurationSpec extends TestKitSuite with FlatSpecLike with Matchers with MockitoSugar with BeforeAndAfter {
  val expectedEndpoint = "oss-cn-shanghai.aliyuncs.com"
  val expectedAccessId = "test-access-id"
  val expectedAccessKey = "test-access-key"
  val expectedToken = Some("test-security-token")
  val expectedFullEndpoint = URI.create("http://oss-cn-shanghai.aliyuncs.com")

  behavior of "TTLOssStorageConfiguration"


  it should "have correct OSS credential info" in {

    val ossConfig = TTLOssStorageConfiguration(TTLOssStorageConfigurationSpec.BcsBackendConfig)

    ossConfig.endpoint shouldEqual expectedEndpoint
    ossConfig.accessId shouldEqual expectedAccessId
    ossConfig.accessKey shouldEqual expectedAccessKey
    ossConfig.securityToken shouldEqual expectedToken

    ossConfig.newOssClient().getEndpoint shouldEqual expectedFullEndpoint

  }
}
