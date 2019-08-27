package cromwell.docker.registryv2.flows.alibabacloudcrregistry

import com.aliyuncs.auth.{BasicCredentials, BasicSessionCredentials}
import cromwell.docker.DockerInfoActor.{DockerInfoContext, DockerInfoFailedResponse, DockerInfoSuccessResponse, DockerInformation}
import cromwell.docker.{DockerHashResult, DockerImageIdentifier, DockerInfoRequest, DockerRegistryConfig}

import net.ceedubs.ficus.Ficus._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.scalatest.mockito.MockitoSugar
import spray.json._

object AlibabaCloudCRRegistrySpec {

  val AlibabaCloudCRRegistryConfigString =
    s"""
       |enable = true
       |# How should docker hashes be looked up. Possible values are "local" and "remote"
       |# "local": Lookup hashes on the local docker daemon using the cli
       |# "remote": Lookup hashes on docker hub and gcr
       |method = "remote"
       |alibabacloudcr {
       |  num-threads = 5
       |  auth {
       |    endpoint = "cr.cn-shanghai.aliyuncs.com"
       |    access-id = "test-access-id"
       |    access-key = "test-access-key"
       |    security-token = "test-security-token"
       |  }
       |}
       |
      """.stripMargin

  val AlibabaCloudCRRegistryConfig = ConfigFactory.parseString(AlibabaCloudCRRegistryConfigString)
}

class AlibabaCloudCRRegistrySpec extends TestKitSuite with FlatSpecLike with Matchers with MockitoSugar with BeforeAndAfter {
  behavior of "AlibabaCloudCRRegistry"

  val hashValue = "fcf39ed78ef0fa27bcc74713b85259alop1b12e6a201e3083af50fd8eda1cbe1"
  val tag = "0.2"
  val notExistTag = "0.3"
  val CRResponse =
    s"""
      |{
      |    "data": {
      |        "total": 2,
      |        "pageSize": 30,
      |        "page": 1,
      |        "tags": [
      |            {
      |                "imageUpdate": 1514432549000,
      |                "imageId": "2842876c9b98f8c7607c1123ks18ff040b76a1d932c6d60c96aa3c283bd221cd",
      |                "digest": "83414d2c3b04e0lo1q7693e31aeca95b82c61949ea8de858579bf16bd92490c6",
      |                "imageSize": 715764,
      |                "tag": "0.1",
      |                "imageCreate": 1514432549000,
      |                "status": "NORMAL"
      |            },
      |            {
      |                "imageUpdate": 1514372113000,
      |                "imageId": "414e6daa772a8cd5dfloqpe503e6e313c372d2e15958ab649709daf9b1065479",
      |                "digest": "$hashValue",
      |                "imageSize": 715653,
      |                "tag": "$tag",
      |                "imageCreate": 1514372044000,
      |                "status": "NORMAL"
      |            }
      |        ]
      |    },
      |    "requestId": "9AFB52D3-6631-4B00-A857-932492097726"
      |}""".stripMargin


  it should "have correct Alibaba Cloud CR image" in {
    val configPath = "alibabacloudcr"
    val registry = new AlibabaCloudCRRegistry(DockerRegistryConfig.fromConfig(AlibabaCloudCRRegistrySpec.AlibabaCloudCRRegistryConfig.as[Config](configPath)).getOrElse(DockerRegistryConfig.default))

    val testCRDockerImage = s"registry.cn-shanghai.aliyuncs.com/batchcompute/ubuntu:$tag"
    val testInvalidCRDockerImage = "registry.cn-not-exist.aliyuncs.com/batchcompute/ubuntu:0.2"
    registry.accepts(DockerImageIdentifier.fromString(testCRDockerImage).get) shouldEqual true
    registry.isValidAlibabaCloudCRHost(Some(testInvalidCRDockerImage)) shouldEqual false
    registry.isValidAlibabaCloudCRHost(None) shouldEqual false
  }

  it should "successfully extract digest from body" in {
    val configPath = "alibabacloudcr"
    val registry = new AlibabaCloudCRRegistry(DockerRegistryConfig.fromConfig(AlibabaCloudCRRegistrySpec.AlibabaCloudCRRegistryConfig.as[Config](configPath)).getOrElse(DockerRegistryConfig.default))

    val testCRDockerImage = s"registry.cn-shanghai.aliyuncs.com/batchcompute/ubuntu:$tag"

    val expectedDockerHashResult = DockerHashResult("sha256", hashValue)
    val expectedDockerInfomation = DockerInformation(expectedDockerHashResult, None)
    val dockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImage).get, List.empty)
    val expectedDockerResponse = DockerInfoSuccessResponse(expectedDockerInfomation, dockerRequest)

    val context: DockerInfoContext = DockerInfoContext(dockerRequest, null)
    registry.extractDigestFromBody(CRResponse.parseJson.asJsObject(), context) shouldEqual expectedDockerResponse
  }

  it should "NOT successfully extract digest from body" in {
    val configPath = "alibabacloudcr"
    val registry = new AlibabaCloudCRRegistry(DockerRegistryConfig.fromConfig(AlibabaCloudCRRegistrySpec.AlibabaCloudCRRegistryConfig.as[Config](configPath)).getOrElse(DockerRegistryConfig.default))

    val testCRDockerImageTagNotExist = s"registry.cn-shanghai.aliyuncs.com/batchcompute/ubuntu:$notExistTag"

    val dockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImageTagNotExist).get, List.empty)
    val context: DockerInfoContext = DockerInfoContext(dockerRequest, null)

    val cRResponseJsObj = CRResponse.parseJson.asJsObject()
    registry.extractDigestFromBody(cRResponseJsObj, context) match {
      case DockerInfoFailedResponse(t, _) => t.getMessage should be(s"Manifest response did not contain a expected tag: $notExistTag, ${cRResponseJsObj}")
      case _ => fail("Failed to get a DockerInfoFailedResponse result.")
    }
  }

  it should "successfully get the correct credentials from context" in {
    val configPath = "alibabacloudcr"
    val registry = new AlibabaCloudCRRegistry(DockerRegistryConfig.fromConfig(AlibabaCloudCRRegistrySpec.AlibabaCloudCRRegistryConfig.as[Config](configPath)).getOrElse(DockerRegistryConfig.default))

    val testCRDockerImageTagNotExist = s"registry.cn-shanghai.aliyuncs.com/batchcompute/ubuntu:$tag"
    val access_id = "test-access-id"
    val access_key = "test-access-key"
    val security_token = "test-token"

    val basicCredential = new BasicCredentials(access_id, access_key)
    val sessionCredential = new BasicSessionCredentials(access_id, access_key, security_token)
    val vpcEndpoint: String = "cr-vpc.cn-shanghai.aliyuncs.com"
    val normalEndpoint = "cr.cn-shanghai.aliyuncs.com"
    val validEndpoint = "cr.validendpoint.com"

    val dockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImageTagNotExist).get, List(basicCredential, normalEndpoint))
    val context: DockerInfoContext = DockerInfoContext(dockerRequest, null)
    registry.getAliyunCredentialFromContext(context) shouldEqual Option(basicCredential)
    registry.getAliyunEndpointFromContext(context) shouldEqual Option(normalEndpoint)

    val vpcDockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImageTagNotExist).get, List(basicCredential, vpcEndpoint))
    val vpcContext: DockerInfoContext = DockerInfoContext(vpcDockerRequest, null)
    registry.getAliyunEndpointFromContext(vpcContext) shouldEqual Option(vpcEndpoint)

    val validDockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImageTagNotExist).get, List(basicCredential, validEndpoint))
    val validContext: DockerInfoContext = DockerInfoContext(validDockerRequest, null)
    registry.getAliyunEndpointFromContext(validContext) shouldEqual None

    val sessionDockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImageTagNotExist).get, List(sessionCredential))
    val sessionContext: DockerInfoContext = DockerInfoContext(sessionDockerRequest, null)
    registry.getAliyunCredentialFromContext(sessionContext) shouldEqual Option(sessionCredential)

    val invalidDockerRequest = DockerInfoRequest(DockerImageIdentifier.fromString(testCRDockerImageTagNotExist).get, List.empty)
    val invalidContext: DockerInfoContext = DockerInfoContext(invalidDockerRequest, null)
    registry.getAliyunCredentialFromContext(invalidContext) shouldEqual None
  }

}
