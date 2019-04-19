package cromwell.backend.google.pipelines.common

import cats.data.NonEmptyList
import cromwell.backend.RuntimeAttributeDefinition
import cromwell.backend.google.pipelines.common.GpuResource.GpuType
import cromwell.backend.google.pipelines.common.PipelinesApiTestConfig.{googleConfiguration, papiAttributes, _}
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk, PipelinesApiWorkingDisk}
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import cromwell.core.WorkflowOptions
import eu.timepit.refined.refineMV
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.helpers.NOPLogger
import org.specs2.mock.Mockito
import spray.json._
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.types._
import wom.values._

class PipelinesApiRuntimeAttributesSpec extends WordSpecLike with Matchers with Mockito {

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]): WorkflowOptions = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  val expectedDefaults = new PipelinesApiRuntimeAttributes(refineMV(1), None, None, Vector("us-central1-b", "us-central1-a"), 0, 10,
    MemorySize(2, MemoryUnit.GB), Vector(PipelinesApiWorkingDisk(DiskType.SSD, 10)), "ubuntu:latest", false,
    ContinueOnReturnCodeSet(Set(0)), false)

  "PipelinesApiRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = Map.empty[String, WomValue]
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "use hardcoded defaults if not declared in task, workflow options, or config (except for docker)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, jesConfiguration = noDefaultsJesConfiguration)
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(failOnStderr = true)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomString("yes"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "validate a valid continueOnReturnCode integer entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomInteger(1))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode boolean entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomBoolean(false))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeFlag(false))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomStringType), Array(WomString("1"), WomString("2"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "validate a valid gpu entry (1)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(1), "gpuType" -> WomString("nvidia-tesla-k80"))
      val expectedRuntimeAttributes = expectedDefaults.copy(gpuResource = Option(GpuResource(gpuCount = refineMV(1), gpuType = GpuType.NVIDIATeslaK80)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid gpu entry (2)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(2), "gpuType" -> WomString("nvidia-tesla-p100"))
      val expectedRuntimeAttributes = expectedDefaults.copy(gpuResource = Option(GpuResource(gpuCount = refineMV(2), gpuType = GpuType.NVIDIATeslaP100)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    // Missing gpu type
    "fail to validate an invalid gpu entry (1)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Please specify a GPU type: nvidia-tesla-p100, nvidia-tesla-k80")
    }

    // Missing gpu count
    "fail to validate an invalid gpu entry (2)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuType" -> WomString("nvidia-tesla-p100"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Please specify how many GPU should be attached to the instance.")
    }

    // unrecoginzed gpu type
    "fail to validate an invalid gpu entry (3)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(1), "gpuType" -> WomString("not-a-gpu"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "not-a-gpu is not a supported GPU type. Supported types are nvidia-tesla-k80, nvidia-tesla-p100")
    }

    // gpu count is not an int
    "fail to validate an invalid gpu entry (4)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting gpuCount runtime attribute to be an Integer")
    }

    "validate a valid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(2))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("2"))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomString("us-central-z"))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-central-z"))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomArray(WomArrayType(WomStringType), Array(WomString("us-central1-y"), WomString("us-central1-z"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-central1-y", "us-central1-z"))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomInteger(3))
      val expectedRuntimeAttributes = expectedDefaults.copy(preemptible = 3)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting preemptible runtime attribute to be an Integer")
    }

    "validate a valid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "bootDiskSizeGb" -> WomInteger(4))
      val expectedRuntimeAttributes = expectedDefaults.copy(bootDiskSize = 4)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "bootDiskSizeGb" -> WomString("4GB"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting bootDiskSizeGb runtime attribute to be an Integer")
    }

    "validate a valid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk 20 SSD"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disks = Seq(PipelinesApiAttachedDisk.parse("local-disk 20 SSD").get))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomInteger(10))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomArray(WomArrayType(WomStringType), Array(WomString("local-disk 20 SSD"), WomString("local-disk 30 SSD"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(disks = Seq(PipelinesApiAttachedDisk.parse("local-disk 20 SSD").get, PipelinesApiAttachedDisk.parse("local-disk 30 SSD").get))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomArray(WomArrayType(WomStringType), Array(WomString("blah"), WomString("blah blah"))))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'")
    }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = MemorySize.parse("1 GB").get)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("blah"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(noAddress = true)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting noAddress runtime attribute to be a Boolean")
    }

    "override config default attributes with default attributes declared in workflow options" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))

      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpu": 2 }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]

      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }

    "override config default runtime attributes with task runtime attributes" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(4))

      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpu": 2 }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]

      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(4))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }

    "override invalid config default attributes with task runtime attributes" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(4))

      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpu": 2.2 }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]

      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(4))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }
  }

  "should parse cpuPlatform correctly" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
    val workflowOptionsJson =
      """{
        |  "default_runtime_attributes": { "cpuPlatform": "the platform" }
        |}
      """.stripMargin.parseJson.asInstanceOf[JsObject]
    val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
    val expectedRuntimeAttributes = expectedDefaults.copy(cpuPlatform = Option("the platform"))
    assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
  }

  private def assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WomValue],
                                                           expectedRuntimeAttributes: PipelinesApiRuntimeAttributes,
                                                           workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                           defaultZones: NonEmptyList[String] = defaultZones,
                                                           jesConfiguration: PipelinesApiConfiguration = papiConfiguration): Unit = {
    try {
      val actualRuntimeAttributes = toJesRuntimeAttributes(runtimeAttributes, workflowOptions, jesConfiguration)
      assert(actualRuntimeAttributes == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertJesRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue],
                                                       exMsg: String,
                                                       workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {
    try {
      toJesRuntimeAttributes(runtimeAttributes, workflowOptions, papiConfiguration)
      fail(s"A RuntimeException was expected with message: $exMsg")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
    ()
  }

  private def toJesRuntimeAttributes(runtimeAttributes: Map[String, WomValue],
                                     workflowOptions: WorkflowOptions,
                                     jesConfiguration: PipelinesApiConfiguration): PipelinesApiRuntimeAttributes = {
    val runtimeAttributesBuilder = PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)
    val defaultedAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(
      staticRuntimeAttributeDefinitions, workflowOptions)(runtimeAttributes)
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)
    PipelinesApiRuntimeAttributes(validatedRuntimeAttributes, jesConfiguration.runtimeConfig)
  }

  private val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get
  private val defaultZones = NonEmptyList.of("us-central1-b", "us-central1-a")
  private val noDefaultsJesConfiguration = new PipelinesApiConfiguration(PipelinesApiTestConfig.NoDefaultsConfigurationDescriptor, genomicsFactory, googleConfiguration, papiAttributes)
  private val staticRuntimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(PipelinesApiTestConfig.papiConfiguration).definitions.toSet
}
