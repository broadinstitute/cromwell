package cromwell.backend.impl.jes

import cromwell.backend.impl.jes.io.{DiskType, JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.{MemorySize, RuntimeAttributeDefinition}
import cromwell.core.WorkflowOptions
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.helpers.NOPLogger
import org.specs2.mock.Mockito
import spray.json._
import wdl4s.parser.MemoryUnit
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString, WdlValue}

class JesRuntimeAttributesSpec extends WordSpecLike with Matchers with Mockito {

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  val expectedDefaults = new JesRuntimeAttributes(1, Vector("us-central1-b"), 0, 10, MemorySize(2, MemoryUnit.GB), Seq(JesWorkingDisk(DiskType.SSD, 10)), None, false, ContinueOnReturnCodeSet(Set(0)), false)
  val expectedDefaultsPlusUbuntuDocker = expectedDefaults.copy(dockerImage = Some("ubuntu:latest"))

  "JesRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = Map.empty[String, WdlValue]
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults.copy(dockerImage = Option("ubuntu:latest"))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WdlInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "failOnStderr" -> WdlBoolean(true))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(failOnStderr = true)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "failOnStderr" -> WdlString("yes"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "validate a valid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlInteger(1))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlArray(WdlArrayType(WdlIntegerType), Array(WdlInteger(1), WdlInteger(2))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlArray(WdlArrayType(WdlStringType), Array(WdlString("1"), WdlString("2"))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "validate a valid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "cpu" -> WdlInteger(2))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(cpu = 2)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "cpu" -> WdlString("2"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(cpu = 2)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "cpu" -> WdlString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid zones entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "zones" -> WdlString("us-central-z"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(zones = Vector("us-central-z"))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid zones entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "zones" -> WdlInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "zones" -> WdlArray(WdlArrayType(WdlStringType), Array(WdlString("us-central1-y"), WdlString("us-central1-z"))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(zones = Vector("us-central1-y", "us-central1-z"))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "zones" -> WdlArray(WdlArrayType(WdlIntegerType), Array(WdlInteger(1), WdlInteger(2))))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "preemptible" -> WdlInteger(3))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(preemptible = 3)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "preemptible" -> WdlString("value"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Failed to validate preemptible runtime attribute: Could not coerce value into an integer")
    }

    "validate a valid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "bootDiskSizeGb" -> WdlInteger(4))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(bootDiskSize = 4)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid bootDiskSizeGb entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "bootDiskSizeGb" -> WdlString("4GB"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting bootDiskSizeGb runtime attribute to be an Integer")
    }

    "validate a valid disks entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "disks" -> WdlString("local-disk 20 SSD"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(disks = Seq(JesAttachedDisk.parse("local-disk 20 SSD").get))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disks entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "disks" -> WdlInteger(10))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "disks" -> WdlArray(WdlArrayType(WdlStringType), Array(WdlString("local-disk 20 SSD"), WdlString("local-disk 30 SSD"))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(disks = Seq(JesAttachedDisk.parse("local-disk 20 SSD").get, JesAttachedDisk.parse("local-disk 30 SSD").get))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "disks" -> WdlArray(WdlArrayType(WdlStringType), Array(WdlString("blah"), WdlString("blah blah"))))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'")
    }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "memory" -> WdlString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(memory = MemorySize.parse("1 GB").get)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "memory" -> WdlString("blah"))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "noAddress" -> WdlBoolean(true))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(noAddress = true)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "noAddress" -> WdlInteger(1))
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Failed to validate noAddress runtime attribute: Could not coerce 1 into a boolean")
    }

    "use reasonable default values" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }
  }

  private def assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue], expectedRuntimeAttributes: JesRuntimeAttributes, workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {
    val withDefaults = RuntimeAttributeDefinition.addDefaultsToAttributes(JesBackendLifecycleActorFactory.staticRuntimeAttributeDefinitions, workflowOptions) _
    try {
      assert(JesRuntimeAttributes(withDefaults(runtimeAttributes), NOPLogger.NOP_LOGGER) == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertJesRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String, workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {
    val withDefaults = RuntimeAttributeDefinition.addDefaultsToAttributes(JesBackendLifecycleActorFactory.staticRuntimeAttributeDefinitions, workflowOptions) _
    try {
      JesRuntimeAttributes(withDefaults(runtimeAttributes), NOPLogger.NOP_LOGGER)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
    ()
  }

  private val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get
}
