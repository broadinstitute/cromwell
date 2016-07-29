package cromwell.backend.sfs

import cromwell.backend.BackendSpec._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation._
import cromwell.core.WorkflowOptions
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.{Logger, LoggerFactory}
import org.specs2.mock.Mockito
import spray.json.{JsArray, JsBoolean, JsNumber, JsObject, JsString, JsValue}
import wdl4s.values.WdlValue

class SharedFileSystemValidatedRuntimeAttributesBuilderSpec extends WordSpecLike with Matchers with Mockito {

  val HelloWorld =
    """
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello ${addressee}!"
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

  val defaultRuntimeAttributes = Map(
    DockerKey -> None,
    FailOnStderrKey -> false,
    ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(0)))

  def workflowOptionsWithDefaultRuntimeAttributes(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map("default_runtime_attributes" -> JsObject(defaults))))
  }

  "SharedFileSystemValidatedRuntimeAttributesBuilder" should {
    "return an instance of itself when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, defaultRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("ubuntu:latest"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "\${addressee}" }""").head
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "use workflow options as default if docker key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("ubuntu:latest"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(Map(DockerKey -> JsString("ubuntu:latest")))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

    "throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      assertRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderrKey -> true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "true" }""").head
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "log a warning and return an instance of itself when tries to validate a valid Docker entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("ubuntu:latest"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      val mockWarnings = new MockWarnings
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        supportsDocker = false, logger = mockWarnings.mockLogger)
      mockWarnings.warnings.size should be(1)
      val message = mockWarnings.warnings.head
      message should include("Unrecognized runtime attribute keys: docker")
    }

    "log a warning and return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (DockerKey -> Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "\${addressee}" }""").head
      val mockWarnings = new MockWarnings
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        supportsDocker = false, logger = mockWarnings.mockLogger)
      mockWarnings.warnings.size should be(1)
      val message = mockWarnings.warnings.head
      message should include("Unrecognized runtime attribute keys: docker")
    }

    "log a warning and throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      val mockWarnings = new MockWarnings
      assertRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String",
        supportsDocker = false, logger = mockWarnings.mockLogger)
      mockWarnings.warnings.size should be(1)
      val message = mockWarnings.warnings.head
      message should include("Unrecognized runtime attribute keys: docker")
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "yes" }""").head
      assertRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderrKey -> true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(Map(FailOnStderrKey -> JsBoolean(true)))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

    "return an instance of itself when tries to validate a valid continueOnReturnCode entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes +
        (ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(1)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: 1 }""").head
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: "value" }""").head
      assertRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting continueOnReturnCode runtime attribute to be either " +
          "a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "use workflow options as default if continueOnReturnCode key is missing" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes +
        (ContinueOnReturnCodeKey -> ContinueOnReturnCodeSet(Set(1, 2)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRuntimeAttributes(
        Map(ContinueOnReturnCodeKey -> JsArray(Vector(JsNumber(1), JsNumber(2)))))
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes,
        workflowOptions = workflowOptions)
    }

    "warn for unrecognized keys" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld,
        """runtime { whatIsThis: "noIdea" andThis: "donno" }""").head
      val mockWarnings = new MockWarnings
      assertRuntimeAttributesSuccessfulCreation(runtimeAttributes, defaultRuntimeAttributes,
        logger = mockWarnings.mockLogger)
      mockWarnings.warnings.size should be(1)
      val message = mockWarnings.warnings.head
      message should include("Unrecognized runtime attribute keys:")
      message should include("whatIsThis")
      message should include("andThis")
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
                                                        supportsDocker: Boolean = true,
                                                        workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                        logger: Logger = defaultLogger): Unit = {

    val builder = SharedFileSystemValidatedRuntimeAttributesBuilder.default.withDockerSupport(supportsDocker)
    val validatedRuntimeAttributes = builder.build(runtimeAttributes, workflowOptions, logger)

    DockerValidation.optional.extract(validatedRuntimeAttributes) should be(
      expectedRuntimeAttributes(DockerKey).asInstanceOf[Option[String]])
    FailOnStderrValidation.default.extract(validatedRuntimeAttributes) should be(
      expectedRuntimeAttributes(FailOnStderrKey).asInstanceOf[Boolean])
    ContinueOnReturnCodeValidation.default.extract(validatedRuntimeAttributes) should be(
      expectedRuntimeAttributes(ContinueOnReturnCodeKey).asInstanceOf[ContinueOnReturnCode])
  }

  private def assertRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String,
                                                    supportsDocker: Boolean = true,
                                                    workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                    logger: Logger = defaultLogger): Unit = {
    val thrown = the[RuntimeException] thrownBy {
      SharedFileSystemValidatedRuntimeAttributesBuilder.default.withDockerSupport(supportsDocker).
        build(runtimeAttributes, workflowOptions, logger)
    }
    thrown.getMessage should include(exMsg)
  }
}
