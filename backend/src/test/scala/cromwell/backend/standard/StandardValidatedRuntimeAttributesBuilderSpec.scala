package cromwell.backend.standard

import cromwell.backend.validation._
import cromwell.backend.{RuntimeAttributeDefinition, TestConfig}
import cromwell.core.WorkflowOptions
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import org.slf4j.{Logger, LoggerFactory}
import org.specs2.mock.Mockito
import spray.json.{JsArray, JsBoolean, JsNumber, JsObject, JsValue}
import wom.RuntimeAttributesKeys._
import wom.values._

class StandardValidatedRuntimeAttributesBuilderSpec extends AnyWordSpecLike with Matchers with Mockito {

  val HelloWorld: String =
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

  def workflowOptionsWithDefaultRuntimeAttributes(defaults: Map[String, JsValue]): WorkflowOptions = {
    WorkflowOptions(JsObject(Map("default_runtime_attributes" -> JsObject(defaults))))
  }

  "SharedFileSystemValidatedRuntimeAttributesBuilder" should {
    "validate when there are no runtime attributes defined" in {
      val runtimeAttributes = Map.empty[String, WomValue]
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, defaultRuntimeAttributes)
    }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("ubuntu:latest"))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      assertRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("failOnStderr" -> WomBoolean(true))
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderrKey -> true)
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "log a warning and validate a valid Docker entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      var warnings = List.empty[Any]
      val mockLogger = mock[Logger]
      mockLogger.warn(anyString).answers(warnings :+= _)
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        includeDockerSupport = false, logger = mockLogger)
      warnings should contain theSameElementsAs List("Unrecognized runtime attribute keys: docker")
    }

    "log a warning and validate an invalid Docker entry" in {
      /*
      NOTE: The behavior used to be: when present, a "docker" runtime attribute would always be validated to ensure that
      the value of the runtime attribute was a String-- even if the actual runtime attribute was unsupported by the
      backend! NOW, if the backend doesn't support docker, and the runtime attributes contain say as an invalid integer,
      say `{ docker: 1 }`, the validation will now succeed. If we want to change this, the git history contains a
      concept of validated-yet-still-unsupported attributes. One could use this behavior to validate other standard
      attributes such as CPU, Memory, etc.-- checking their syntax even when they are unsupported by the current
      backend.

      https://github.com/broadinstitute/cromwell/blob/a4a952de33f6d1ef646be51d298c3d613a8cce5f/supportedBackends/sfs/src/main/scala/cromwell/backend/sfs/SharedFileSystemValidatedRuntimeAttributesBuilder.scala
       */
      val expectedRuntimeAttributes = defaultRuntimeAttributes
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      var warnings = List.empty[Any]
      val mockLogger = mock[Logger]
      mockLogger.warn(anyString).answers(warnings :+= _)
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        includeDockerSupport = false, logger = mockLogger)
      warnings should contain theSameElementsAs List("Unrecognized runtime attribute keys: docker")
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("failOnStderr" -> WomString("yes"))
      assertRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderrKey -> true)
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(Map(FailOnStderrKey -> JsBoolean(true)))
      val runtimeAttributes = Map.empty[String, WomValue]
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

    "validate a valid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("continueOnReturnCode" -> WomInteger(1))
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(1)))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("continueOnReturnCode" -> WomString("value"))
      assertRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "use workflow options as default if continueOnReturnCode key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes +
        (ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(1, 2)))
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(
        Map(ContinueOnReturnCodeKey -> JsArray(Vector(JsNumber(1), JsNumber(2)))))
      val runtimeAttributes = Map.empty[String, WomValue]
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

  }

  val defaultLogger: Logger = LoggerFactory.getLogger(classOf[StandardValidatedRuntimeAttributesBuilderSpec])
  val emptyWorkflowOptions: WorkflowOptions = WorkflowOptions.fromMap(Map.empty).get

  val mockBackendRuntimeConfig = Option(TestConfig.optionalRuntimeConfig)

  private def assertRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WomValue],
                                                        expectedRuntimeAttributes: Map[String, Any],
                                                        includeDockerSupport: Boolean = true,
                                                        workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                        logger: Logger = defaultLogger): Unit = {

    val builder = if (includeDockerSupport) {
      StandardValidatedRuntimeAttributesBuilder.default(mockBackendRuntimeConfig).withValidation(DockerValidation.optional)
    } else {
      StandardValidatedRuntimeAttributesBuilder.default(mockBackendRuntimeConfig)
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

  private def assertRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue], exMsg: String,
                                                    supportsDocker: Boolean = true,
                                                    workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                    logger: Logger = defaultLogger): Unit = {
    val thrown = the[RuntimeException] thrownBy {
      val builder = if (supportsDocker) {
        StandardValidatedRuntimeAttributesBuilder.default(mockBackendRuntimeConfig).withValidation(DockerValidation.optional)
      } else {
        StandardValidatedRuntimeAttributesBuilder.default(mockBackendRuntimeConfig)
      }
      val runtimeAttributeDefinitions = builder.definitions.toSet
      val addDefaultsToAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributeDefinitions, workflowOptions) _

      builder.build(addDefaultsToAttributes(runtimeAttributes), logger)
    }
    thrown.getMessage should include(exMsg)
    ()
  }
}
