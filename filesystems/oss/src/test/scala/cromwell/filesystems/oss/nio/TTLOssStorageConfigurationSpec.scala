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

//  before {
//    BcsMount.pathBuilders = List(mockPathBuiler)
//  }
//
//  val jobId = "test-bcs-job"
//  val mockOssConf = DefaultOssStorageConfiguration("oss.aliyuncs.com", "test-id", "test-key")
//  val mockPathBuiler = OssPathBuilder(mockOssConf)
//  val mockPathBuilders = List(mockPathBuiler)
//  lazy val workflowDescriptor =  buildWdlWorkflowDescriptor(
//    SampleWdl.HelloWorld.workflowSource(),
//    inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
//  )
//  lazy val jobKey = {
//    val call = workflowDescriptor.callable.taskCallNodes.head
//    BackendJobDescriptorKey(call, None, 1)
//  }
//
//
//  val expectedContinueOnReturn = ContinueOnReturnCodeSet(Set(0))
//  val expectedDockerTag = Some(BcsDockerWithPath("ubuntu/latest", "oss://bcs-reg/ubuntu/"))
//  val expectedDocker = Some(BcsDockerWithoutPath("registry.cn-beijing.aliyuncs.com/test/testubuntu:0.1"))
//  val expectedFailOnStderr = false
//  val expectedUserData = Some(Vector(new BcsUserData("key", "value")))
//  //val expectedMounts = Some(Vector(new BcsInputMount(Left(mockPathBuiler.build("oss://bcs-bucket/bcs-dir/").get), Left(DefaultPathBuilder.build("/home/inputs/").get), false)))
//  val expectedMounts = Some(Vector(new BcsInputMount(Left(mockPathBuiler.build("oss://bcs-bucket/bcs-dir/").get), Right("/home/inputs/"), false)))
//  val expectedCluster = Some(Left("cls-mycluster"))
//  val expectedSystemDisk = Some(BcsSystemDisk("cloud", 50))
//  val expectedDataDsik = Some(BcsDataDisk("cloud", 250, "/home/data/"))
//
//  val expectedReserveOnFail = Some(true)
//  val expectedAutoRelease = Some(true)
//  val expectedWorkerPath = Some("oss://bcs-bucket/workflow/worker.tar.gz")
//  val expectedTimeout = Some(3000)
//  val expectedVerbose = Some(false)
//  val expectedVpc = Some(BcsVpcConfiguration(Some("192.168.0.0/16"), Some("vpc-xxxx")))
//  val expectedTag = Some("jobTag")
//
//
//  val expectedRuntimeAttributes = new BcsRuntimeAttributes(expectedContinueOnReturn, expectedDockerTag, expectedDocker, expectedFailOnStderr,  expectedMounts, expectedUserData, expectedCluster,
//    expectedSystemDisk, expectedDataDsik, expectedReserveOnFail, expectedAutoRelease, expectedWorkerPath, expectedTimeout, expectedVerbose, expectedVpc, expectedTag)
//
//
//  protected def createBcsRuntimeAttributes(runtimeAttributes: Map[String, WomValue]): BcsRuntimeAttributes = {
//    val builder = BcsRuntimeAttributes.runtimeAttributesBuilder(BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
//    val default = RuntimeAttributeDefinition.addDefaultsToAttributes(
//      builder.definitions.toSet, BcsTestUtilSpec.EmptyWorkflowOption)(runtimeAttributes)
//    val validated = builder.build(default, NOPLogger.NOP_LOGGER)
//    BcsRuntimeAttributes(validated, BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
//  }

  val expectedEndpoint = "oss-cn-shanghai.aliyuncs.com"
  val expectedAccessId = "test-access-id"
  val expectedAccessKey = "test-access-key"
  val expectedToken = Some("test-security-token")
  val expectedFullEndpoint = URI.create("http://oss-cn-shanghai.aliyuncs.com")

  behavior of s"TTLOssStorageConfiguration"


  it should "have correct OSS credential info" in {

    val ossConfig = TTLOssStorageConfiguration(TTLOssStorageConfigurationSpec.BcsBackendConfig)

    ossConfig.endpoint shouldEqual expectedEndpoint
    ossConfig.accessId shouldEqual expectedAccessId
    ossConfig.accessKey shouldEqual expectedAccessKey
    ossConfig.securityToken shouldEqual expectedToken

    ossConfig.newOssClient().getEndpoint shouldEqual expectedFullEndpoint

  }
}
