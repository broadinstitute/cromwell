package cromwell.backend.impl.tes

import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.{BackendConfigurationDescriptor, RuntimeAttributeDefinition, TestConfig}
import cromwell.core.WorkflowOptions
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.helpers.NOPLogger
import spray.json._
import squants.information.Gigabytes
import wom.types._
import wom.values._

class TesRuntimeAttributesSpec extends WordSpecLike with Matchers {

  val expectedDefaults = new TesRuntimeAttributes(
    ContinueOnReturnCodeSet(Set(0)),
    "ubuntu:latest",
    None,
    false,
    None,
    None,
    None
  )

  val expectedDefaultsPlusUbuntuDocker = expectedDefaults.copy(dockerImage = "ubuntu:latest")

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

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
      assertFailure(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "validate a valid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomInteger(1))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomStringType), Array(WomString("1"), WomString("2"))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomString("value"))
      assertFailure(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
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
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = Option(Gigabytes(1)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("blah"))
      assertFailure(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid disk entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disk" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disk = Option(Gigabytes(1)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disk entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disk" -> WomString("blah"))
      assertFailure(runtimeAttributes, "Expecting disk runtime attribute to be an Integer or String with format '8 GB'")
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
  }

  private val mockConfigurationDescriptor = BackendConfigurationDescriptor(TesTestConfig.backendConfig, TestConfig.globalConfig)
  private val mockTesConfiguration = new TesConfiguration(mockConfigurationDescriptor)

  private def assertSuccess(runtimeAttributes: Map[String, WomValue],
                            expectedRuntimeAttributes: TesRuntimeAttributes,
                            workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {

    try {
      val actualRuntimeAttributes = toTesRuntimeAttributes(runtimeAttributes, workflowOptions, mockTesConfiguration)
      assert(actualRuntimeAttributes == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertFailure(runtimeAttributes: Map[String, WomValue],
                            exMsg: String,
                            workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {
    try {
      toTesRuntimeAttributes(runtimeAttributes, workflowOptions, mockTesConfiguration)
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
                                     tesConfiguration: TesConfiguration): TesRuntimeAttributes = {
    val runtimeAttributesBuilder = TesRuntimeAttributes.runtimeAttributesBuilder(tesConfiguration.runtimeConfig)
    val defaultedAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(
      staticRuntimeAttributeDefinitions, workflowOptions)(runtimeAttributes)
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)
    TesRuntimeAttributes(validatedRuntimeAttributes, tesConfiguration.runtimeConfig
    )
  }
}
