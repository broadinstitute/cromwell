package cromwell.backend.impl.bcs

import com.typesafe.config.ConfigFactory
import common.collections.EnhancedCollections._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, RuntimeAttributeDefinition}
import cromwell.backend.BackendSpec.buildWdlWorkflowDescriptor
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.oss.OssPathBuilder
import cromwell.filesystems.oss.nio.DefaultOssStorageConfiguration
import cromwell.util.SampleWdl
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.scalatest.mockito.MockitoSugar
import org.slf4j.helpers.NOPLogger
import spray.json.{JsObject, JsString}
import wom.values.WomValue

object BcsTestUtilSpec {

  val DefaultRunAttributesString =
    """
      |default-runtime-attributes {
      |  failOnStderr: false
      |  continueOnReturnCode: 0
      |  cluster: "cls-mycluster"
      |  mounts: "oss://bcs-bucket/bcs-dir/ /home/inputs/ false"
      |  dockerTag: "ubuntu/latest oss://bcs-reg/ubuntu/"
      |  docker: "registry.cn-beijing.aliyuncs.com/test/testubuntu:0.1"
      |  userData: "key value"
      |  reserveOnFail: true
      |  autoReleaseJob: true
      |  verbose: false
      |  systemDisk: "cloud 50"
      |  dataDisk: "cloud 250 /home/data/"
      |  timeout: 3000
      |  vpc: "192.168.0.0/16 vpc-xxxx"
      |  tag: "jobTag"
      |  imageId: "img-ubuntu-vpc"
      |}
    """.stripMargin

  val BcsBackendConfigString =
    s"""
      |root = "oss://your-bucket/cromwell-exe"
      |dockerRoot = "/cromwell-executions"
      |region = ""
      |
      |access-id = ""
      |access-key = ""
      |security-token = ""
      |
      |filesystems {
      |  oss {
      |    auth {
      |        endpoint = ""
      |        access-id = ""
      |        access-key = ""
      |        security-token = ""
      |    }
      |    caching {
      |        duplication-strategy = "reference"
      |    }
      |  }
      |}
      |
      |$DefaultRunAttributesString
      |
      |""".stripMargin

  val BcsBackendConfigWithoutDefaultString =
    s"""
       |root = "oss://your-bucket/cromwell-exe"
       |dockerRoot = "/cromwell-executions"
       |region = ""
       |
       |access-id = ""
       |access-key = ""
       |security-token = ""
       |
       |filesystems {
       |  oss {
       |    auth {
       |        endpoint = ""
       |        access-id = ""
       |        access-key = ""
       |        security-token = ""
       |    }
       |  }
       |}
       |
       |""".stripMargin

  val BcsGlobalConfigString =
    s"""
      |backend {
      |  default = "BCS"
      |  providers {
      |    BCS {
      |      actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
      |      config {
      |      $BcsBackendConfigString
      |      }
      |    }
      |  }
      |}
      |
      |""".stripMargin

  val BcsBackendConfig = ConfigFactory.parseString(BcsBackendConfigString)
  val BcsGlobalConfig = ConfigFactory.parseString(BcsGlobalConfigString)
  val BcsBackendConfigWithoutDefault = ConfigFactory.parseString(BcsBackendConfigWithoutDefaultString)
  val BcsBackendConfigurationDescriptor = BackendConfigurationDescriptor(BcsBackendConfig, BcsGlobalConfig)
  val BcsBackendConfigurationWithoutDefaultDescriptor = BackendConfigurationDescriptor(BcsBackendConfigWithoutDefault, BcsGlobalConfig)
  val EmptyWorkflowOption = WorkflowOptions.fromMap(Map.empty).get
}

trait BcsTestUtilSpec extends TestKitSuite with FlatSpecLike with Matchers with MockitoSugar with BeforeAndAfter {

  before {
    BcsMount.pathBuilders = List(mockPathBuilder)
  }

  val jobId = "test-bcs-job"
  val mockOssConf = DefaultOssStorageConfiguration("oss.aliyuncs.com", "test-id", "test-key")
  val mockPathBuilder = OssPathBuilder(mockOssConf)
  val mockPathBuilders = List(mockPathBuilder)
  lazy val workflowDescriptor =  buildWdlWorkflowDescriptor(
    SampleWdl.HelloWorld.workflowSource(),
    inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
  )
  lazy val jobKey = {
    val call = workflowDescriptor.callable.taskCallNodes.head
    BackendJobDescriptorKey(call, None, 1)
  }


  val expectedContinueOnReturn = ContinueOnReturnCodeSet(Set(0))
  val expectedDockerTag = Option(BcsDockerWithPath("ubuntu/latest", "oss://bcs-reg/ubuntu/"))
  val expectedDocker = Option(BcsDockerWithoutPath("registry.cn-beijing.aliyuncs.com/test/testubuntu:0.1"))
  val expectedFailOnStderr = false
  val expectedUserData = Option(Vector(new BcsUserData("key", "value")))
  val expectedMounts = Option(Vector(new BcsInputMount(Left(mockPathBuilder.build("oss://bcs-bucket/bcs-dir/").get), Right("/home/inputs/"), false)))
  val expectedCluster = Option(Left("cls-mycluster"))
  val expectedImageId = Option("img-ubuntu-vpc")
  val expectedSystemDisk = Option(BcsSystemDisk("cloud", 50))
  val expectedDataDisk = Option(BcsDataDisk("cloud", 250, "/home/data/"))

  val expectedReserveOnFail = Option(true)
  val expectedAutoRelease = Option(true)
  val expectedTimeout = Option(3000)
  val expectedVerbose = Option(false)
  val expectedVpc = Option(BcsVpcConfiguration(Option("192.168.0.0/16"), Option("vpc-xxxx")))
  val expectedTag = Option("jobTag")


  val expectedRuntimeAttributes = new BcsRuntimeAttributes(expectedContinueOnReturn, expectedDockerTag, expectedDocker, expectedFailOnStderr,  expectedMounts, expectedUserData, expectedCluster,
    expectedImageId, expectedSystemDisk, expectedDataDisk, expectedReserveOnFail, expectedAutoRelease, expectedTimeout, expectedVerbose, expectedVpc, expectedTag)


  protected def createBcsRuntimeAttributes(runtimeAttributes: Map[String, WomValue]): BcsRuntimeAttributes = {
    val builder = BcsRuntimeAttributes.runtimeAttributesBuilder(BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
    val default = RuntimeAttributeDefinition.addDefaultsToAttributes(
      builder.definitions.toSet, BcsTestUtilSpec.EmptyWorkflowOption)(runtimeAttributes)
    val validated = builder.build(default, NOPLogger.NOP_LOGGER)
    BcsRuntimeAttributes(validated, BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
  }
}
