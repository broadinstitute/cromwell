package cromwell.backend.google.batch.models

import cats.data.NonEmptyList
import cromwell.backend.RuntimeAttributeDefinition
import cromwell.backend.google.batch.models.GcpBatchTestConfig._
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.google.batch.io.{DiskType, GcpBatchAttachedDisk, GcpBatchWorkingDisk}
import cromwell.core.WorkflowOptions
import eu.timepit.refined.refineMV
import org.scalatest.TestSuite
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import org.slf4j.helpers.NOPLogger
import spray.json._
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}

final class GcpBatchRuntimeAttributesSpec
    extends AnyWordSpecLike
    with Matchers
    with GcpBatchRuntimeAttributesSpecsMixin {

  "GcpBatchRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = Map.empty[String, WomValue]
      assertBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "use hardcoded defaults if not declared in task, workflow options, or config (except for docker)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes,
                                                     expectedRuntimeAttributes,
                                                     batchConfiguration = noDefaultsBatchConfiguration
      )
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      assertBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(failOnStderr = true)
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomString("yes"))
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'"
      )
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomString("value"))
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"
      )
    }

    "validate a valid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(2))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("2"))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV(2))
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("value"))
      assertBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomString("us-central-z"))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-central-z"))
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomInteger(1))
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]"
      )
    }

    "validate a valid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"),
                                  "zones" -> WomArray(WomArrayType(WomStringType),
                                                      List(WomString("us-central1-y"), WomString("us-central1-z"))
                                  )
      )
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-central1-y", "us-central1-z"))
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"),
                                  "zones" -> WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2)))
      )
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]"
      )
    }

    "validate a valid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomInteger(3))
      val expectedRuntimeAttributes = expectedDefaults.copy(preemptible = 3)
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomString("value"))
      assertBatchRuntimeAttributesFailedCreation(runtimeAttributes,
                                                 "Expecting preemptible runtime attribute to be an Integer"
      )
    }

    "validate a valid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "bootDiskSizeGb" -> WomInteger(4))
      val expectedRuntimeAttributes = expectedDefaults.copy(bootDiskSize = 4)
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "bootDiskSizeGb" -> WomString("4GB"))
      assertBatchRuntimeAttributesFailedCreation(runtimeAttributes,
                                                 "Expecting bootDiskSizeGb runtime attribute to be an Integer"
      )
    }

    "validate a valid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk 20 SSD"))
      val expectedRuntimeAttributes =
        expectedDefaults.copy(disks = Seq(GcpBatchAttachedDisk.parse("local-disk 20 SSD").get))
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomInteger(10))
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Expecting disks runtime attribute to be a comma separated String or Array[String]"
      )
    }

    "fail to validate a valid disks array entry" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"),
            "disks" -> WomArray(WomArrayType(WomStringType), List(WomString("blah"), WomString("blah blah")))
        )
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'"
      )
    }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = MemorySize.parse("1 GB").get)
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("blah"))
      assertBatchRuntimeAttributesFailedCreation(
        runtimeAttributes,
        "Expecting memory runtime attribute to be an Integer or String with format '8 GB'"
      )
    }

    "validate a valid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(noAddress = true)
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomInteger(1))
      assertBatchRuntimeAttributesFailedCreation(runtimeAttributes,
                                                 "Expecting noAddress runtime attribute to be a Boolean"
      )
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
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
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
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
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
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
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
      assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }
  }
}

trait GcpBatchRuntimeAttributesSpecsMixin {
  this: TestSuite =>

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]): WorkflowOptions =
    WorkflowOptions(
      JsObject(
        Map(
          "default_runtime_attributes" -> JsObject(defaults)
        )
      )
    )

  val expectedDefaults = new GcpBatchRuntimeAttributes(
    cpu = refineMV(1),
    cpuPlatform = None,
    gpuResource = None,
    zones = Vector("us-central1-b", "us-central1-a"),
    preemptible = 0,
    bootDiskSize = 10,
    memory = MemorySize(2, MemoryUnit.GB),
    disks = Vector(GcpBatchWorkingDisk(DiskType.SSD, 10)),
    dockerImage = "ubuntu:latest",
    failOnStderr = false,
    continueOnReturnCode = ContinueOnReturnCodeSet(Set(0)),
    noAddress = false,
    useDockerImageCache = None,
    checkpointFilename = None
  )

  def assertBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WomValue],
                                                     expectedRuntimeAttributes: GcpBatchRuntimeAttributes,
                                                     workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                     defaultZones: NonEmptyList[String] = defaultZones,
                                                     batchConfiguration: GcpBatchConfiguration = gcpBatchConfiguration
  ): Unit = {
    try {
      val actualRuntimeAttributes = toBatchRuntimeAttributes(runtimeAttributes, workflowOptions, batchConfiguration)
      assert(actualRuntimeAttributes == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  def assertBatchRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue],
                                                 exMsgs: List[String],
                                                 workflowOptions: WorkflowOptions
  ): Unit = {
    Try(toBatchRuntimeAttributes(runtimeAttributes, workflowOptions, gcpBatchConfiguration)) match {
      case Success(oops) =>
        fail(
          s"Expected error containing strings: ${exMsgs.map(s => s"'$s'").mkString(", ")} but instead got Success($oops)"
        )
      case Failure(ex) => exMsgs foreach { exMsg => assert(ex.getMessage.contains(exMsg)) }
    }
    ()
  }

  def assertBatchRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue],
                                                 exMsg: String,
                                                 workflowOptions: WorkflowOptions = emptyWorkflowOptions
  ): Unit =
    assertBatchRuntimeAttributesFailedCreation(runtimeAttributes, List(exMsg), workflowOptions)

  def toBatchRuntimeAttributes(runtimeAttributes: Map[String, WomValue],
                               workflowOptions: WorkflowOptions,
                               batchConfiguration: GcpBatchConfiguration
  ): GcpBatchRuntimeAttributes = {
    val runtimeAttributesBuilder = GcpBatchRuntimeAttributes.runtimeAttributesBuilder(batchConfiguration)
    val defaultedAttributes =
      RuntimeAttributeDefinition.addDefaultsToAttributes(staticRuntimeAttributeDefinitions, workflowOptions)(
        runtimeAttributes
      )
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)
    GcpBatchRuntimeAttributes(validatedRuntimeAttributes, batchConfiguration.runtimeConfig)
  }

  val emptyWorkflowOptions: WorkflowOptions = WorkflowOptions.fromMap(Map.empty).get
  val defaultZones: NonEmptyList[String] = NonEmptyList.of("us-central1-b", "us-central1-a")
  val noDefaultsBatchConfiguration = new GcpBatchConfiguration(GcpBatchTestConfig.NoDefaultsConfigurationDescriptor,
                                                               googleConfiguration,
                                                               batchAttributes
  )
  val staticRuntimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    GcpBatchRuntimeAttributes.runtimeAttributesBuilder(GcpBatchTestConfig.gcpBatchConfiguration).definitions.toSet
}
