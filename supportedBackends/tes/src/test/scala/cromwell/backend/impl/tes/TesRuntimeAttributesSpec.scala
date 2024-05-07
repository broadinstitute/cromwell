package cromwell.backend.impl.tes

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.{BackendConfigurationDescriptor, RuntimeAttributeDefinition, TestConfig}
import cromwell.core.WorkflowOptions
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import org.slf4j.helpers.NOPLogger
import spray.json._
import wom.format.MemorySize
import wom.types._
import wom.values._

class TesRuntimeAttributesSpec extends AnyWordSpecLike with CromwellTimeoutSpec with Matchers {

  val expectedDefaults = new TesRuntimeAttributes(
    ContinueOnReturnCodeSet(Set(0)),
    "ubuntu:latest",
    None,
    false,
    None,
    None,
    None,
    false,
    None,
    Map.empty
  )

  val expectedDefaultsPlusUbuntuDocker = expectedDefaults.copy(dockerImage = "ubuntu:latest")

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) =
    WorkflowOptions(
      JsObject(
        Map(
          "default_runtime_attributes" -> JsObject(defaults)
        )
      )
    )

  "TesRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = Map.empty[String, WomValue]
      assertFailure(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults.copy(dockerImage = "ubuntu:latest")
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      assertFailure(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(failOnStderr = true)
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomString("yes"))
      assertFailure(
        runtimeAttributes,
        "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'"
      )
    }

    "validate a valid preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(preemptible = true)
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid azureSasEnvironmentVariable entry" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"), TesRuntimeAttributes.LocalizedSasKey -> WomString("THIS_IS_VALID"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(localizedSasEnvVar = Some("THIS_IS_VALID"))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid azureSasEnvironmentVariable entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"),
                                  TesRuntimeAttributes.LocalizedSasKey -> WomString("THIS IS INVALID")
      )
      assertFailure(runtimeAttributes, "Value must be a string containing only letters, numbers, and underscores.")
    }

    "convert a positive integer preemptible entry to true boolean" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomInteger(3))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(preemptible = true)
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "convert a nonpositive integer preemptible entry to false boolean" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomInteger(-1))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(preemptible = false)
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "convert a valid string preemptible entry to boolean" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomString("true"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(preemptible = true)
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid string preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomString("yes"))
      assertFailure(
        runtimeAttributes,
        "Expecting preemptible runtime attribute to be an Integer, Boolean, or a String with values of 'true' or 'false'"
      )
    }

    "fail to validate an invalid type preemptible entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "preemptible" -> WomFloat(3.14))
      assertFailure(
        runtimeAttributes,
        "Expecting preemptible runtime attribute to be an Integer, Boolean, or a String with values of 'true' or 'false'"
      )
    }

    "validate a valid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomInteger(1))
      val expectedRuntimeAttributes =
        expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"),
            "continueOnReturnCode" -> WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2)))
        )
      val expectedRuntimeAttributes =
        expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"),
            "continueOnReturnCode" -> WomArray(WomArrayType(WomStringType), List(WomString("1"), WomString("2")))
        )
      val expectedRuntimeAttributes =
        expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomString("value"))
      assertFailure(
        runtimeAttributes,
        "Expecting returnCodes/continueOnReturnCode" +
          " runtime attribute to be either a String '*', 'true', or 'false', a Boolean, or an Array[Int]."
      )
    }

    "validate a valid cpu entry" in assertSuccess(
      Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(2)),
      expectedDefaultsPlusUbuntuDocker.copy(cpu = Option(refineMV[Positive](2)))
    )

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("2"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(cpu = Option(refineMV[Positive](2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("value"))
      assertFailure(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = Option(MemorySize.parse("1 GB").get))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("blah"))
      assertFailure(runtimeAttributes,
                    "Expecting memory runtime attribute to be an Integer or String with format '8 GB'"
      )
    }

    "validate a valid disk entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disk" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disk = Option(MemorySize.parse("1 GB").get))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disk entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disk" -> WomString("blah"))
      assertFailure(runtimeAttributes, "Expecting disk runtime attribute to be an Integer or String with format '8 GB'")
    }

    "parse an HDD definition" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk 10 HDD"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disk = Option(MemorySize.parse("10 GB").get))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "parse an SSD definition" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk 10 SSD"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disk = Option(MemorySize.parse("10 GB").get))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "refuse multiple `local-disk` instances" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk 10 SSD, local-disk 20 SSD"))
      assertFailure(runtimeAttributes, "Expecting exactly one disk definition on this backend, found multiple")
    }

    "refuse custom mount points" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("/some/mnt 20 SSD"))
      assertFailure(runtimeAttributes, "Disks with custom mount points are not supported by this backend")
    }

    "refuse custom AND multiple mount points" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("/mnt/tmp 10 LOCAL, local-disk 20 HDD"))
      assertFailure(runtimeAttributes, "Disks with custom mount points are not supported by this backend")
    }

    "not accept a single comma" ignore {
      // Surprisingly, the PAPI code we call under the covers validates `,` and give the user a default disk.
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString(","))
      assertFailure(
        runtimeAttributes,
        "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got: ','"
      )
    }

    "not accept empty string" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString(""))
      assertFailure(
        runtimeAttributes,
        "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got: ''"
      )
    }

    "not accept `banana`" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("banana"))
      assertFailure(
        runtimeAttributes,
        "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got: 'banana'"
      )
    }

    "not accept a random number (chosen by fair dice roll)" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomInteger(4))
      assertFailure(runtimeAttributes,
                    "Expecting disks runtime attribute to be a comma separated String or Array[String]"
      )
    }

    "validate a valid dockerWorkingDir entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "dockerWorkingDir" -> WomString("/tmp"))
      val expectedRuntimeAttributes = expectedDefaults.copy(dockerWorkingDir = Option("/tmp"))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid dockerWorkingDir entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "dockerWorkingDir" -> WomInteger(1))
      assertFailure(runtimeAttributes, "Expecting dockerWorkingDir runtime attribute to be a String")
    }

    "use reasonable default values" in assertSuccess(
      Map("docker" -> WomString("ubuntu:latest")),
      expectedDefaultsPlusUbuntuDocker
    )

    "not turn unknown string attributes into backend parameters when using default config" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "foo" -> WomString("bar"))
      assertSuccess(runtimeAttributes, expectedDefaults)
    }

    "turn unknown string attributes into backend parameters" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "foo" -> WomString("bar"))
      val expectedRuntimeAttributes = expectedDefaults.copy(backendParameters = Map("foo" -> Option("bar")))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes, tesConfig = mockTesConfigWithBackendParams)
    }

    "exclude unknown non-string attributes from backend parameters" in {
      val runtimeAttributes =
        Map("docker" -> WomString("ubuntu:latest"), "foo" -> WomInteger(5), "bar" -> WomString("baz"))
      val expectedRuntimeAttributes = expectedDefaults.copy(backendParameters = Map("bar" -> Option("baz")))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes, tesConfig = mockTesConfigWithBackendParams)
    }

    "turn populated optional unknown string attributes into backend parameters" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "foo" -> WomOptionalValue(WomString("bar")))
      val expectedRuntimeAttributes = expectedDefaults.copy(backendParameters = Map("foo" -> Option("bar")))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes, tesConfig = mockTesConfigWithBackendParams)
    }

    "turn unpopulated optional unknown string attributes into backend parameters" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "foo" -> WomOptionalValue.none(WomStringType))
      val expectedRuntimeAttributes = expectedDefaults.copy(backendParameters = Map("foo" -> None))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes, tesConfig = mockTesConfigWithBackendParams)
    }
  }

  private val mockConfigurationDescriptor =
    BackendConfigurationDescriptor(TesTestConfig.backendConfig, TestConfig.globalConfig)
  private val mockTesConfiguration = new TesConfiguration(mockConfigurationDescriptor)
  private val mockTesConfigWithBackendParams = new TesConfiguration(
    mockConfigurationDescriptor.copy(backendConfig = TesTestConfig.backendConfigWithBackendParams)
  )

  private def assertSuccess(runtimeAttributes: Map[String, WomValue],
                            expectedRuntimeAttributes: TesRuntimeAttributes,
                            workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                            tesConfig: TesConfiguration = mockTesConfiguration
  ): Unit = {

    try {
      val actualRuntimeAttributes = toTesRuntimeAttributes(runtimeAttributes, workflowOptions, tesConfig)
      assert(actualRuntimeAttributes == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertFailure(runtimeAttributes: Map[String, WomValue],
                            exMsg: String,
                            workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                            tesConfig: TesConfiguration = mockTesConfiguration
  ): Unit = {
    try {
      toTesRuntimeAttributes(runtimeAttributes, workflowOptions, tesConfig)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
    ()
  }

  private val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get
  private val staticRuntimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    TesRuntimeAttributes.runtimeAttributesBuilder(mockTesConfiguration.runtimeConfig).definitions.toSet

  private def toTesRuntimeAttributes(runtimeAttributes: Map[String, WomValue],
                                     workflowOptions: WorkflowOptions,
                                     tesConfiguration: TesConfiguration
  ): TesRuntimeAttributes = {
    val runtimeAttributesBuilder = TesRuntimeAttributes.runtimeAttributesBuilder(tesConfiguration.runtimeConfig)
    val defaultedAttributes =
      RuntimeAttributeDefinition.addDefaultsToAttributes(staticRuntimeAttributeDefinitions, workflowOptions)(
        runtimeAttributes
      )
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)
    TesRuntimeAttributes(validatedRuntimeAttributes, runtimeAttributes, tesConfiguration)
  }
}
