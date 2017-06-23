package cromwell.backend.impl.tes

import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.{BackendConfigurationDescriptor, MemorySize, RuntimeAttributeDefinition, TestConfig}
import cromwell.core.WorkflowOptions
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.helpers.NOPLogger
import spray.json._
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString, WdlValue}

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
      val runtimeAttributes = Map.empty[String, WdlValue]
      assertFailure(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults.copy(dockerImage = "ubuntu:latest")
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WdlInteger(1))
      assertFailure(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "failOnStderr" -> WdlBoolean(true))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(failOnStderr = true)
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "failOnStderr" -> WdlString("yes"))
      assertFailure(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "validate a valid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlInteger(1))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlArray(WdlArrayType(WdlIntegerType), Array(WdlInteger(1), WdlInteger(2))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlArray(WdlArrayType(WdlStringType), Array(WdlString("1"), WdlString("2"))))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "continueOnReturnCode" -> WdlString("value"))
      assertFailure(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "validate a valid cpu entry" in assertSuccess(
      Map("docker" -> WdlString("ubuntu:latest"), "cpu" -> WdlInteger(2)),
      expectedDefaultsPlusUbuntuDocker.copy(cpu = Option(2))
    )

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "cpu" -> WdlString("2"))
      val expectedRuntimeAttributes = expectedDefaultsPlusUbuntuDocker.copy(cpu = Option(2))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "cpu" -> WdlString("value"))
      assertFailure(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "memory" -> WdlString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = Option(MemorySize.parse("1 GB").get))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "memory" -> WdlString("blah"))
      assertFailure(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid disk entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "disk" -> WdlString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disk = Option(MemorySize.parse("1 GB").get))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disk entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "disk" -> WdlString("blah"))
      assertFailure(runtimeAttributes, "Expecting disk runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid dockerWorkingDir entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "dockerWorkingDir" -> WdlString("/tmp"))
      val expectedRuntimeAttributes = expectedDefaults.copy(dockerWorkingDir = Option("/tmp"))
      assertSuccess(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid dockerWorkingDir entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"), "dockerWorkingDir" -> WdlInteger(1))
      assertFailure(runtimeAttributes, "Expecting dockerWorkingDir runtime attribute to be a String")
    }

    "use reasonable default values" in assertSuccess(
      Map("docker" -> WdlString("ubuntu:latest")),
      expectedDefaultsPlusUbuntuDocker
    )
  }

  private val mockConfigurationDescriptor = BackendConfigurationDescriptor(TesTestConfig.backendConfig, TestConfig.globalConfig)
  private val mockTesConfiguration = new TesConfiguration(mockConfigurationDescriptor)

  private def assertSuccess(runtimeAttributes: Map[String, WdlValue],
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

  private def assertFailure(runtimeAttributes: Map[String, WdlValue],
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


  private def toTesRuntimeAttributes(runtimeAttributes: Map[String, WdlValue],
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
