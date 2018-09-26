package cromwell.backend.impl.bcs

import com.typesafe.config.ConfigFactory
import common.collections.EnhancedCollections._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, RuntimeAttributeDefinition}
import cromwell.backend.BackendSpec.{buildWdlWorkflowDescriptor}
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.oss.OssPathBuilder
import cromwell.filesystems.oss.nio.OssStorageConfiguration
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
      |  docker: "ubuntu/latest oss://bcs-reg/ubuntu/"
      |  userData: "key value"
      |  reserveOnFail: true
      |  autoReleaseJob: true
      |  verbose: false
      |  workerPath: "oss://bcs-bucket/workflow/worker.tar.gz"
      |  systemDisk: "cloud 50"
      |  dataDisk: "cloud 250 /home/data/"
      |  timeout: 3000
      |  vpc: "192.168.0.0/16 vpc-xxxx"
      |  tag: "jobTag"
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
    BcsMount.pathBuilders = List(mockPathBuiler)
  }

  val jobId = "test-bcs-job"
  val mockOssConf = OssStorageConfiguration("oss.aliyuncs.com", "test-id", "test-key")
  val mockPathBuiler = OssPathBuilder(mockOssConf)
  val mockPathBuilders = List(mockPathBuiler)
  lazy val workflowDescriptor =  buildWdlWorkflowDescriptor(
    SampleWdl.HelloWorld.workflowSource(),
    inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
  )
  lazy val jobKey = {
    val call = workflowDescriptor.callable.taskCallNodes.head
    BackendJobDescriptorKey(call, None, 1)
  }


  val expectedContinueOnReturn = ContinueOnReturnCodeSet(Set(0))
  val expectedDocker = Some(BcsDockerWithPath("ubuntu/latest", "oss://bcs-reg/ubuntu/"))
  val expectedFailOnStderr = false
  val expectedUserData = Some(Vector(new BcsUserData("key", "value")))
  val expectedMounts = Some(Vector(new BcsInputMount(mockPathBuiler.build("oss://bcs-bucket/bcs-dir/").get, DefaultPathBuilder.build("/home/inputs/").get, false)))
  val expectedCluster = Some(Left("cls-mycluster"))
  val expectedSystemDisk = Some(BcsSystemDisk("cloud", 50))
  val expectedDataDsik = Some(BcsDataDisk("cloud", 250, "/home/data/"))

  val expectedReserveOnFail = Some(true)
  val expectedAutoRelease = Some(true)
  val expectedWorkerPath = Some("oss://bcs-bucket/workflow/worker.tar.gz")
  val expectedTimeout = Some(3000)
  val expectedVerbose = Some(false)
  val expectedVpc = Some(BcsVpcConfiguration(Some("192.168.0.0/16"), Some("vpc-xxxx")))
  val expectedTag = Some("jobTag")


  val expectedRuntimeAttributes = new BcsRuntimeAttributes(expectedContinueOnReturn, expectedDocker, expectedFailOnStderr,  expectedMounts, expectedUserData, expectedCluster,
    expectedSystemDisk, expectedDataDsik, expectedReserveOnFail, expectedAutoRelease, expectedWorkerPath, expectedTimeout, expectedVerbose, expectedVpc, expectedTag)


  protected def createBcsRuntimeAttributes(runtimeAttributes: Map[String, WomValue]): BcsRuntimeAttributes = {
    val builder = BcsRuntimeAttributes.runtimeAttributesBuilder(BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
    val default = RuntimeAttributeDefinition.addDefaultsToAttributes(
      builder.definitions.toSet, BcsTestUtilSpec.EmptyWorkflowOption)(runtimeAttributes)
    val validated = builder.build(default, NOPLogger.NOP_LOGGER)
    BcsRuntimeAttributes(validated, BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
  }
}
