package cromwell.engine.workflow

import akka.actor.ActorRef
import akka.testkit.TestDuration
import cromwell.CromwellTestkitSpec
import cromwell.core.{WorkflowOptions, WorkflowId}
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.{Backend, BackendType, CromwellBackend}
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor.{MaterializationFailure, MaterializationSuccess, MaterializeWorkflow}
import cromwell.util.SampleWdl.HelloWorld
import org.mockito.Mockito._
import org.scalatest.BeforeAndAfter
import org.scalatest.mock.MockitoSugar
import spray.json.DefaultJsonProtocol._
import spray.json._
import wdl4s.ThrowableWithErrors

import scala.concurrent.duration._
import scala.language.postfixOps

class MaterializeWorkflowDescriptorActorSpec
  extends CromwellTestkitSpec with BeforeAndAfter with MockitoSugar {
  val MissingInputsJson = "{}"
  val MalformedJson = "foobar bad json!"
  val MalformedWdl = "foobar bad wdl!"

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
  val Timeout = 1.second.dilated
  var materializeWfActor: ActorRef = _
  var validWorkflowSources = WorkflowSourceFiles(
    validWdlSource, validInputsJson, customizedLocalBackendOptions)

  before {
    materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props())
    // Needed since we might want to run this test as test-only
    CromwellBackend.initBackends(List("local"), "local", system)
    CromwellBackend.registerCustomBackend("retryableCallsSpecBackend", backendMock)
  }

  after {
    system.stop(materializeWfActor)
    reset(backendMock)
  }

  "MaterializeWorkflowDescriptorActor" should {
    "return MaterializationSuccess when Namespace, options, runtime attributes and inputs are valid" in {
      within(Timeout) {
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(), validWorkflowSources)
        expectMsgPF() {
          case MaterializationSuccess(wfDesc) =>
            if (wfDesc.namespace.tasks.size != 1) fail("Number of tasks is not equals to one.")
            if (!wfDesc.coercedInputs.head._1.contains("hello.hello.addressee")) fail("Input does not contains 'hello.hello.addressee'.")
            if (!wfDesc.workflowOptions.get("backend").get.contains("retryableCallsSpecBackend"))
              fail("Workflow option does not comply with 'backend' entry.")
          case unknown =>
            fail(s"Response is not equal to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return MaterializationFailure when there is an invalid Namespace coming from a WDL source" in {
      within(Timeout) {
        val malformedSources =
          validWorkflowSources.copy(wdlSource = MalformedWdl)
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(), malformedSources)
        expectMsgPF() {
          case MaterializationFailure(failure) =>
            failure match {
              case validationException: IllegalArgumentException with ThrowableWithErrors =>
                if (!validationException.message.contains("Workflow input processing failed."))
                  fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Unable to load namespace from workflow"))
                  fail("Message from error nr 1 in validation exception does not contains 'Unable to load namespace from workflow'.")
            }
          case unknown => fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return MaterializationFailure when there are workflow options coming from a malformed workflow options JSON file" in {
      within(Timeout) {
        val malformedSources =
          validWorkflowSources.copy(workflowOptionsJson = malformedOptionsJson)
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(), malformedSources)
        expectMsgPF() {
          case MaterializationFailure(failure) =>
            failure match {
              case validationException: IllegalArgumentException with ThrowableWithErrors =>
                if (!validationException.message.contains("Workflow input processing failed."))
                  fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Workflow contains invalid options JSON"))
                  fail("Message from error nr 1 in validation exception does not contains 'Workflow contains invalid options JSON'.")
            }
          case unknown =>
            fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(0)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return MaterializationFailure when there are invalid workflow options coming from a workflow options JSON file" in {
      when(backendMock.assertWorkflowOptions(
        WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get))
        .thenThrow(new IllegalStateException("Some exception"))
      within(Timeout) {
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(), validWorkflowSources)
        expectMsgPF() {
          case MaterializationFailure(failure) =>
            failure match {
              case validationException: IllegalArgumentException with ThrowableWithErrors =>
                if (!validationException.message.contains("Workflow input processing failed."))
                  fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1)  fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Workflow has invalid options for backend"))
                  fail("Message from error nr 1 in validation exception does not contains 'Workflow has invalid options for backend'.")
            }
          case unknown =>
            fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return MaterializationFailure when there are workflow inputs coming from a malformed workflow inputs JSON file" in {
      within(Timeout) {
        val malformedSources =
          validWorkflowSources.copy(inputsJson = malformedInputsJson)
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(), malformedSources)
        expectMsgPF() {
          case MaterializationFailure(failure) =>
            failure match {
              case validationException: IllegalArgumentException with ThrowableWithErrors =>
                if (!validationException.message.contains("Workflow input processing failed."))
                  fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Workflow contains invalid inputs JSON"))
                  fail("Message from error nr 1 in validation exception does not contains 'Workflow contains invalid inputs JSON'")
            }
          case unknown =>
            fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return MaterializationFailure when there are invalid workflow inputs coming from a workflow inputs JSON file" in {
      within(Timeout) {
        val malformedSources =
          validWorkflowSources.copy(inputsJson = invalidInputJson)
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(),
          malformedSources)
        expectMsgPF() {
          case MaterializationFailure(failure) =>
            failure match {
              case validationException: IllegalArgumentException with ThrowableWithErrors =>
                if (!validationException.message.contains("Workflow input processing failed."))
                  fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1) fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Required workflow input 'hello.hello.addressee' not specified."))
                  fail("Message from error nr 1 in validation exception does not contains 'Required workflow input 'hello.hello.addressee' not specified'.")
            }
          case unknown =>
            fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(
        WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }

    "return MaterializationFailure when there are invalid runtime requirements coming from WDL file" in {
      when(backendMock.backendType).thenReturn(BackendType.JES)
      within(Timeout) {
        val malformedSources =
          validWorkflowSources.copy(wdlSource = invalidWdlSource)
        materializeWfActor ! MaterializeWorkflow(WorkflowId.randomId(), malformedSources)
        expectMsgPF() {
          case MaterializationFailure(failure) =>
            failure match {
              case validationException: IllegalArgumentException with ThrowableWithErrors =>
                if (!validationException.message.contains("Workflow input processing failed."))
                  fail("Message coming from validation exception does not contains 'Workflow input processing failed'.")
                if (validationException.errors.size != 1)  fail("Number of errors coming from validation exception is not equals to one.")
                if (!validationException.errors.head.contains("Failed to validate runtime attributes."))
                  fail("Message from error nr 1 in validation exception does not contains 'Failed to validate runtime attributes'.")
            }
          case unknown =>
            fail(s"Response is not equals to the expected one. Response: $unknown")
        }
      }
      verify(backendMock, times(1)).assertWorkflowOptions(WorkflowOptions.fromJsonString(customizedLocalBackendOptions).get)
    }
  }
}
