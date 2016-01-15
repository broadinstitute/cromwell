package cromwell

import akka.event.LoggingAdapter
import akka.testkit.{EventFilter, TestActorRef}
import cromwell.engine._
import cromwell.engine.backend.BackendType
import cromwell.engine.backend.jes.{JesBackend, JesBackendCall}
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.workflow.{CallKey, WorkflowManagerActor}
import cromwell.logging.WorkflowLogger
import cromwell.util.SampleWdl
import org.slf4j.LoggerFactory
import org.specs2.mock.Mockito
import wdl4s.{RuntimeAttributes, ThrowableWithErrors, CallInputs}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Try, Success}

class InvalidRuntimeAttributesSpec extends CromwellTestkitSpec with Mockito {

  "A workflow with a task with invalid runtime attributes" should {
    "fail on JES Backend" in {
      val wd = mock[WorkflowDescriptor]
      val logger = LoggerFactory.getLogger("Fake logger")
      wd.workflowLogger returns logger
      wd.shortId returns "fakeShortID"

      val jesBackendCall = mock[JesBackendCall]
      def exception = {
        Try(CromwellRuntimeAttributes(RuntimeAttributes(Map.empty[String, Seq[String]]), BackendType.JES)).failed.get
      }

      jesBackendCall.runtimeAttributes throws exception
      jesBackendCall.workflowLoggerWithCall(any[String], any[Option[LoggingAdapter]]) answers { (params, mock) =>
        val array: Array[Object] = params.asInstanceOf[Array[Object]]
        WorkflowLogger("WorkflowActor", wd, akkaLogger = array(1).asInstanceOf[Option[LoggingAdapter]])
      }

      val mockJes = mock[JesBackend]
      mockJes.initializeForWorkflow(any[WorkflowDescriptor]) returns Success(Map.empty)
      mockJes.cleanUpForWorkflow(any[WorkflowDescriptor])(any[ExecutionContext]) returns Future.successful(())
      mockJes.bindCall(any[WorkflowDescriptor], any[CallKey], any[CallInputs],any[AbortRegistrationFunction]) returns jesBackendCall

      val workflowSources = WorkflowSourceFiles(SampleWdl.HelloWorld.wdlSource(), SampleWdl.HelloWorld.wdlJson, "{}")
      val submitMessage = WorkflowManagerActor.SubmitWorkflow(workflowSources)
      runWdlWithWorkflowManagerActor(
        wma = TestActorRef(new WorkflowManagerActor(mockJes)),
        submitMsg = submitMessage,
        stdout = Map.empty,
        stderr = Map.empty,
        eventFilter = EventFilter.error(pattern = "RuntimeAttribute is not valid.", occurrences = 1),
        terminalState = WorkflowFailed
      )
    }

    "fail on Local Backend" in {
      runWdl(
        sampleWdl = SampleWdl.HelloWorld,
        runtime =
          """ runtime { wrongAttribute: "nop" }""".stripMargin,
        eventFilter = EventFilter.error(pattern = ".*RuntimeAttribute is not valid.*", occurrences = 1),
        terminalState = WorkflowFailed
      )
    }
  }

}
