package cromwell.webservice

import java.util.UUID

import cromwell.util.SampleWdl.HelloWorld

import akka.actor.{Actor, Props}
import cromwell.binding._
import cromwell.engine.WorkflowManagerActor.{WorkflowStatus, SubmitWorkflow}
import cromwell.engine._
import org.scalatest.{Matchers, FlatSpec}
import spray.http._
import spray.testkit.ScalatestRouteTest

import spray.json._
import DefaultJsonProtocol._

object MockWorkflowManagerActor {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowRawInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage

  val createdWorkflowId = UUID.randomUUID()
  val runningWorkflowId = UUID.randomUUID()
  val unknownId = UUID.randomUUID()

  val submittedWorkflowId = UUID.randomUUID()

  def props: Props = Props(classOf[MockWorkflowManagerActor])
}

class MockWorkflowManagerActor extends Actor  {

  def receive = {
    case SubmitWorkflow(wdl, inputs) =>
      sender ! MockWorkflowManagerActor.submittedWorkflowId

    case WorkflowStatus(id) =>
      val msg = id match {
        case MockWorkflowManagerActor.runningWorkflowId =>
          Some(WorkflowRunning)
        case _ =>
          None
      }
      sender ! msg
  }
}

object CromwellApiServiceSpec {
  val MalformedInputsJson : String = "foobar bad json!"
  val MalformedWdl : String = "foobar bad wdl!"
}

class CromwellApiServiceSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {
  def actorRefFactory = system
  val workflowManagerActorRef = system.actorOf(Props(new MockWorkflowManagerActor()))

  "CromwellApiService" should "return 404 for get of unknown workflow" in {
    Get(s"/workflow/${MockWorkflowManagerActor.unknownId}") ~>
      sealRoute(queryRoute) ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for get of a malformed workflow id" in {
    Get(s"/workflow/foobar/status") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return 200 for get of a known workflow id" in {
    Get(s"/workflow/${MockWorkflowManagerActor.runningWorkflowId}/status") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }

        assertResult( s"""{"id":"${MockWorkflowManagerActor.runningWorkflowId.toString}","status":"WorkflowRunning"}""") {
          responseAs[String]
        }
      }
  }


  "Cromwell submit workflow API" should "return 201 for a succesful workfow submission " in {
    Post("/workflows", FormData(Seq("wdlSource" -> HelloWorld.WdlSource, "workflowInputs" -> HelloWorld.RawInputs.toJson.toString()))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.Created) {
          status
        }
        assertResult( s"""{"id":"${MockWorkflowManagerActor.submittedWorkflowId.toString}","status":"WorkflowSubmitted"}""") {
          responseAs[String]
        }
      }
  }

  it should "return 400 for a malformed JSON " in {
    Post("/workflows", FormData(Seq("wdlSource" -> HelloWorld.WdlSource, "workflowInputs" -> CromwellApiServiceSpec.MalformedInputsJson))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult("workflowInput JSON was malformed") {
          responseAs[String]
        }
      }
  }

}
