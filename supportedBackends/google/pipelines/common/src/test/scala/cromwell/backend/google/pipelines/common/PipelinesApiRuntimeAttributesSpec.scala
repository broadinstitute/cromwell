package cromwell.backend.google.pipelines.common

import cats.data.NonEmptyList
import cromwell.backend.RuntimeAttributeDefinition
import cromwell.backend.google.pipelines.common.PipelinesApiTestConfig.{googleConfiguration, papiAttributes, _}
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk, PipelinesApiWorkingDisk}
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import cromwell.core.WorkflowOptions
import eu.timepit.refined.refineMV
import org.scalatest.{Matchers, TestSuite, WordSpecLike}
import org.slf4j.helpers.NOPLogger
import org.specs2.mock.Mockito
import spray.json._
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}

final class PipelinesApiRuntimeAttributesSpec
  extends WordSpecLike
    with Matchers
    with Mockito
    with PipelinesApiRuntimeAttributesSpecsMixin {

  "PipelinesApiRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = Map.empty[String, WomValue]
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "use hardcoded defaults if not declared in task, workflow options, or config (except for docker)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, papiConfiguration = noDefaultsPapiConfiguration)
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(failOnStderr = true)
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomString("yes"))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "validate a valid continueOnReturnCode integer entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomInteger(1))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode boolean entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomBoolean(false))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeFlag(false))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomStringType), Array(WomString("1"), WomString("2"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomString("value"))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "validate a valid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(2))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("2"))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("value"))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomString("us-central-z"))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-central-z"))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomInteger(1))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomArray(WomArrayType(WomStringType), Array(WomString("us-central1-y"), WomString("us-central1-z"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-central1-y", "us-central1-z"))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomInteger(3))
      val expectedRuntimeAttributes = expectedDefaults.copy(preemptible = 3)
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomString("value"))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting preemptible runtime attribute to be an Integer")
    }

    "validate a valid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "bootDiskSizeGb" -> WomInteger(4))
      val expectedRuntimeAttributes = expectedDefaults.copy(bootDiskSize = 4)
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "bootDiskSizeGb" -> WomString("4GB"))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting bootDiskSizeGb runtime attribute to be an Integer")
    }

    "validate a valid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk 20 SSD"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disks = Seq(PipelinesApiAttachedDisk.parse("local-disk 20 SSD").get))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomInteger(10))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomArray(WomArrayType(WomStringType), Array(WomString("local-disk 20 SSD"), WomString("local-disk 30 SSD"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(disks = Seq(PipelinesApiAttachedDisk.parse("local-disk 20 SSD").get, PipelinesApiAttachedDisk.parse("local-disk 30 SSD").get))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomArray(WomArrayType(WomStringType), Array(WomString("blah"), WomString("blah blah"))))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'")
    }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = MemorySize.parse("1 GB").get)
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("blah"))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(noAddress = true)
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomInteger(1))
      assertPapiRuntimeAttributesFailedCreation(runtimeAttributes,
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
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
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
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
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
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }

    "parse cpuPlatform correctly" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpuPlatform": "the platform" }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]
      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpuPlatform = Option("the platform"))
      assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }
  }
}

trait PipelinesApiRuntimeAttributesSpecsMixin { this: TestSuite =>

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]): WorkflowOptions = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  val expectedDefaults = new PipelinesApiRuntimeAttributes(refineMV(1), None, None, Vector("us-central1-b", "us-central1-a"), 0, 10,
    MemorySize(2, MemoryUnit.GB), Vector(PipelinesApiWorkingDisk(DiskType.SSD, 10)), "ubuntu:latest", false,
    ContinueOnReturnCodeSet(Set(0)), false)

  def assertPapiRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WomValue],
                                                    expectedRuntimeAttributes: PipelinesApiRuntimeAttributes,
                                                    workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                    defaultZones: NonEmptyList[String] = defaultZones,
                                                    papiConfiguration: PipelinesApiConfiguration = papiConfiguration): Unit = {
    try {
      val actualRuntimeAttributes = toPapiRuntimeAttributes(runtimeAttributes, workflowOptions, papiConfiguration)
      assert(actualRuntimeAttributes == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  def assertPapiRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue],
                                                exMsgs: List[String],
                                                workflowOptions: WorkflowOptions): Unit = {
    Try(toPapiRuntimeAttributes(runtimeAttributes, workflowOptions, papiConfiguration)) match {
      case Success(oops) =>
        fail(s"Expected error containing strings: ${exMsgs.map(s => s"'$s'").mkString(", ")} but instead got Success($oops)")
      case Failure(ex) => exMsgs foreach { exMsg =>  assert(ex.getMessage.contains(exMsg)) }
    }
    ()
  }

  def assertPapiRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue],
                                                exMsg: String,
                                                workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {
    assertPapiRuntimeAttributesFailedCreation(runtimeAttributes, List(exMsg), workflowOptions)
  }

  def toPapiRuntimeAttributes(runtimeAttributes: Map[String, WomValue],
                              workflowOptions: WorkflowOptions,
                              papiConfiguration: PipelinesApiConfiguration): PipelinesApiRuntimeAttributes = {
    val runtimeAttributesBuilder = PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(papiConfiguration)
    val defaultedAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(
      staticRuntimeAttributeDefinitions, workflowOptions)(runtimeAttributes)
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)
    PipelinesApiRuntimeAttributes(validatedRuntimeAttributes, papiConfiguration.runtimeConfig)
  }

  val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get
  val defaultZones = NonEmptyList.of("us-central1-b", "us-central1-a")
  val noDefaultsPapiConfiguration = new PipelinesApiConfiguration(PipelinesApiTestConfig.NoDefaultsConfigurationDescriptor, genomicsFactory, googleConfiguration, papiAttributes)
  val staticRuntimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(PipelinesApiTestConfig.papiConfiguration).definitions.toSet
}
