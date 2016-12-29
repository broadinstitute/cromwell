package cromwell.backend.sfs

import cromwell.backend.RuntimeAttributeDefinition
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation._
import cromwell.core.WorkflowOptions
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.{Logger, LoggerFactory}
import org.specs2.mock.Mockito
import spray.json.{JsArray, JsBoolean, JsNumber, JsObject, JsValue}
import wdl4s.values.{WdlBoolean, WdlInteger, WdlString, WdlValue}

class SharedFileSystemValidatedRuntimeAttributesBuilderSpec extends WordSpecLike with Matchers with Mockito {

  val HelloWorld =
    s"""
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello $${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin


  val defaultRuntimeAttributes: Map[String, Any] = Map(
    DockerKey -> None,
    FailOnStderrKey -> false,
    ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(0)))

  def workflowOptionsWithDefaultRuntimeAttributes(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map("default_runtime_attributes" -> JsObject(defaults))))
  }

  "SharedFileSystemValidatedRuntimeAttributesBuilder" should {
    "validate when there are no runtime attributes defined" in {
      val runtimeAttributes = Map.empty[String, WdlValue]
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, defaultRuntimeAttributes)
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"))
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("ubuntu:latest"))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WdlInteger(1))
      assertRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("failOnStderr" -> WdlBoolean(true))
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderrKey -> true)
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "log a warning and validate a valid Docker entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes
      val runtimeAttributes = Map("docker" -> WdlString("ubuntu:latest"))
      val mockWarnings = new MockWarnings
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        includeDockerSupport = false, logger = mockWarnings.mockLogger)
      mockWarnings.warnings.size should be(1)
      val message = mockWarnings.warnings.head
      message should include("Unrecognized runtime attribute keys: docker")
    }

    "log a warning and fail to validate an invalid Docker entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes
      val runtimeAttributes = Map("docker" -> WdlInteger(1))
      val mockWarnings = new MockWarnings
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        includeDockerSupport = false, logger = mockWarnings.mockLogger)
      mockWarnings.warnings.size should be(1)
      val message = mockWarnings.warnings.head
      message should include("Unrecognized runtime attribute keys: docker")
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("failOnStderr" -> WdlString("yes"))
      assertRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderrKey -> true)
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(Map(FailOnStderrKey -> JsBoolean(true)))
      val runtimeAttributes = Map.empty[String, WdlValue]
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

    "validate a valid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("continueOnReturnCode" -> WdlInteger(1))
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(1)))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("continueOnReturnCode" -> WdlString("value"))
      assertRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "use workflow options as default if continueOnReturnCode key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes +
        (ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(1, 2)))
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(
        Map(ContinueOnReturnCodeKey -> JsArray(Vector(JsNumber(1), JsNumber(2)))))
      val runtimeAttributes = Map.empty[String, WdlValue]
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

    /**
      * Captures warning log messages with a mock.
      *
      * TODO: A generic MockLogger (with all levels captured) may be an Slf4J+specs2 replacement for the Slf4J+Logback
      * TestLogger in lenthall.
      */
    class MockWarnings() {
      var warnings: Seq[String] = Seq.empty
      val mockLogger = mock[Logger]
      mockLogger.warn(anyString).answers { result =>
        result match {
          case message: String =>
            warnings :+= message
        }
      }
    }
  }

  val defaultLogger = LoggerFactory.getLogger(classOf[SharedFileSystemValidatedRuntimeAttributesBuilderSpec])
  val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get

  private def assertRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue],
                                                        expectedRuntimeAttributes: Map[String, Any],
                                                        includeDockerSupport: Boolean = true,
                                                        workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                        logger: Logger = defaultLogger): Unit = {

    val builder = if (includeDockerSupport) {
      SharedFileSystemValidatedRuntimeAttributesBuilder.default.withValidation(DockerValidation.optional)
    } else {
      SharedFileSystemValidatedRuntimeAttributesBuilder.default
    }
    val runtimeAttributeDefinitions = builder.definitions.toSet
    val addDefaultsToAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, workflowOptions) _

    val validatedRuntimeAttributes = builder.build(addDefaultsToAttributes(runtimeAttributes), logger)

    val docker = RuntimeAttributesValidation.extractOption(
      DockerValidation.instance, validatedRuntimeAttributes)
    val failOnStderr = RuntimeAttributesValidation.extract(
      FailOnStderrValidation.instance, validatedRuntimeAttributes)
    val continueOnReturnCode = RuntimeAttributesValidation.extract(
      ContinueOnReturnCodeValidation.instance, validatedRuntimeAttributes)

    docker should be(expectedRuntimeAttributes(DockerKey).asInstanceOf[Option[String]])
    failOnStderr should be(expectedRuntimeAttributes(FailOnStderrKey).asInstanceOf[Boolean])
    continueOnReturnCode should be(
      expectedRuntimeAttributes(ContinueOnReturnCodeKey).asInstanceOf[ContinueOnReturnCode])
    ()
  }

  private def assertRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String,
                                                    supportsDocker: Boolean = true,
                                                    workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                    logger: Logger = defaultLogger): Unit = {
    val thrown = the[RuntimeException] thrownBy {
      val builder = if (supportsDocker) {
        SharedFileSystemValidatedRuntimeAttributesBuilder.default.withValidation(DockerValidation.optional)
      } else {
        SharedFileSystemValidatedRuntimeAttributesBuilder.default
      }
      val runtimeAttributeDefinitions = builder.definitions.toSet
      val addDefaultsToAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, workflowOptions) _

      builder.build(addDefaultsToAttributes(runtimeAttributes), logger)
    }
    thrown.getMessage should include(exMsg)
    ()
  }
}
