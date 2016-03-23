package cromwell.engine.workflow

import akka.actor.{ActorRef, ActorSystem}
import akka.testkit.{DefaultTimeout, ImplicitSender, TestDuration, TestKit}
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.{Backend, BackendType, CromwellBackend}
import cromwell.engine.workflow.ValidateActor.{ValidateWorkflow, ValidationFailure, ValidationSuccess}
import cromwell.util.SampleWdl.HelloWorld
import org.mockito.Mockito._
import org.scalatest.mock.MockitoSugar
import org.scalatest.{BeforeAndAfter, Matchers, WordSpecLike}
import spray.json.DefaultJsonProtocol._
import spray.json._
import wdl4s.ValidationException

import scala.concurrent.duration._
import scala.language.postfixOps

class ValidateActorSpec extends TestKit(ActorSystem("ValidateActorSpecSystem"))
  with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfter with MockitoSugar {
  val MissingInputsJson: String = "{}"
  val MalformedJson: String = "foobar bad json!"
  val MalformedWdl: String = "foobar bad wdl!"

  val customizedLocalBackendOptions =
    """
      |{
      |  "backend": "retryableCallsSpecBackend"
      |}""".stripMargin

  val malformedOptionsJson = MalformedJson
  val validInputsJson = HelloWorld.rawInputs.toJson.toString()
  val malformedInputsJson = MalformedJson
  val invalidInputJson = Map("llo.addsee" -> "wod").toJson.toString()
  val validRunAttr =
    """
      |runtime {
      |  docker: "ubuntu:latest"
      |}
    """.stripMargin
  val invalidRunAttr =
    """
      |runtime {
      |  asda: "sda"
      |}
    """.stripMargin
  val validWdlSource = HelloWorld.wdlSource(validRunAttr)
  val invalidWdlSource = HelloWorld.wdlSource(invalidRunAttr)

  val backendMock = mock[Backend]

  var validateActor: ActorRef = _

  val Timeout = 1.second.dilated

  before {
    validateActor = system.actorOf(ValidateActor.props())
    // Needed since we might want to run this test as test-only
    CromwellBackend.initBackends(List("local"), "local", system)
    CromwellBackend.registerCustomBackend("retryableCallsSpecBackend", backendMock)
  }

  after {
    system.stop(validateActor)
    reset(backendMock)
  }

  "A Validate actor" should {
    "return ValidationSuccess when Namespace, options, runtime attributes and inputs are valid" in {
      within(Timeout) {
        validateActor ! ValidateWorkflow(validWdlSource, Some(validInputsJson), Some(customizedLocalBackendOptions))
        expectMsgPF() {
          case validationSuccess: ValidationSuccess =>
            if (validationSuccess.namespaceWithWorkflow.tasks.size != 1) fail("Number of tasks is not equals to one.")
            if (!validationSuccess.coercedInputs.get.head._1.contains("hello.hello.addressee")) fail("Input does not contains 'hello.hello.addressee'.")
            if (!validationSuccess.workflowOptions.get.get("backend").get.contains("retryableCallsSpecBackend")) fail("Workflow option does not comply with 'backend' entry.")
            if (!validationSuccess.runtimeAttributes.head.head.contains("docker")) fail("Runtime attribute does not contains 'docker' entry.")
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return ValidationFailure when there is an invalid Namespace coming from a WDL source" in {
      within(Timeout) {
        validateActor ! ValidateWorkflow(MalformedWdl, Some(validInputsJson), Some(customizedLocalBackendOptions))
        expectMsgPF() {
          case validationFailure: ValidationFailure =>
            validationFailure.reason match {
              case validationException: ValidationException =>
                if (!validationException.message.contains("Workflow input processing failed.")) fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Unable to load namespace from workflow")) fail("Message from error nr 1 in validation exception does not contains 'Unable to load namespace from workflow'.")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return ValidationFailure when there are workflow options coming from a malformed workflow options JSON file" in {
      within(Timeout) {
        validateActor ! ValidateWorkflow(validWdlSource, Some(validInputsJson), Some(malformedOptionsJson))
        expectMsgPF() {
          case validationFailure: ValidationFailure =>
            validationFailure.reason match {
              case validationException: ValidationException =>
                if (!validationException.message.contains("Workflow input processing failed.")) fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Workflow contains invalid options JSON")) fail("Message from error nr 1 in validation exception does not contains 'Workflow contains invalid options JSON'.")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(0)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return ValidationFailure when there are invalid workflow options coming from a workflow options JSON file" in {
      when(backendMock.assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)).thenThrow(new IllegalStateException("Some exception"))
      within(Timeout) {
        validateActor ! ValidateWorkflow(validWdlSource, Some(validInputsJson), Some(customizedLocalBackendOptions))
        expectMsgPF() {
          case validationFailure: ValidationFailure =>
            validationFailure.reason match {
              case validationException: ValidationException =>
                if (!validationException.message.contains("Workflow input processing failed.")) fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Workflow has invalid options for backend")) fail("Message from error nr 1 in validation exception does not contains 'Workflow has invalid options for backend'.")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return ValidationFailure when there are workflow inputs coming from a malformed workflow inputs JSON file" in {
      within(Timeout) {
        validateActor ! ValidateWorkflow(validWdlSource, Some(malformedInputsJson), Some(customizedLocalBackendOptions))
        expectMsgPF() {
          case validationFailure: ValidationFailure =>
            validationFailure.reason match {
              case validationException: ValidationException =>
                if (!validationException.message.contains("Workflow input processing failed.")) fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Workflow contains invalid inputs JSON")) fail("Message from error nr 1 in validation exception does not contains 'Workflow contains invalid inputs JSON'")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return ValidationFailure when there are invalid workflow inputs coming from a workflow inputs JSON file" in {
      within(Timeout) {
        validateActor ! ValidateWorkflow(validWdlSource, Some(invalidInputJson), Some(customizedLocalBackendOptions))
        expectMsgPF() {
          case validationFailure: ValidationFailure =>
            validationFailure.reason match {
              case validationException: ValidationException =>
                if (!validationException.message.contains("Workflow input processing failed.")) fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Required workflow input 'hello.hello.addressee' not specified.")) fail("Message from error nr 1 in validation exception does not contains 'Required workflow input 'hello.hello.addressee' not specified'.")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return ValidationFailure when there are invalid runtime requirements coming from WDL file" in {
      when(backendMock.backendType).thenReturn(BackendType.JES)
      within(Timeout) {
        validateActor ! ValidateWorkflow(invalidWdlSource, Some(validInputsJson), Some(customizedLocalBackendOptions))
        expectMsgPF() {
          case validationFailure: ValidationFailure =>
            validationFailure.reason match {
              case validationException: ValidationException =>
                if (!validationException.message.contains("Workflow input processing failed.")) fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Failed to validate runtime attributes.")) fail("Message from error nr 1 in validation exception does not contains 'Failed to validate runtime attributes'.")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }
  }
}
