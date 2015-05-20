package cromwell.webservice

import java.util.UUID

import scala.concurrent.ExecutionContext.Implicits.global
import akka.actor.{Actor, Props}
import cromwell.binding._
import cromwell.engine.WorkflowManagerActor.{WorkflowStatus, SubmitWorkflow}
import cromwell.engine._
import org.scalatest.{Matchers, FlatSpec}
import spray.http._
import spray.testkit.ScalatestRouteTest
import akka.pattern.pipe

import scala.concurrent.Future

object MockWorkflowManagerActor {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage

  val runningWorkflowId = UUID.randomUUID()
  val unknownId = UUID.randomUUID()

  def props: Props = Props(classOf[MockWorkflowManagerActor])
}

class MockWorkflowManagerActor extends Actor  {

  def receive = {
    case SubmitWorkflow(wdl, inputs) =>
      println("not yet implemented")

    case WorkflowStatus(id) =>
      Future {
        id match {
          case MockWorkflowManagerActor.runningWorkflowId =>
            Some(WorkflowRunning)
          case _ =>
            None
        }
      }.pipeTo(sender())

  }
}

class CromwellApiServiceSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {

  def actorRefFactory = system
  val workflowManagerActorRef = system.actorOf(Props(new MockWorkflowManagerActor()))
  val unknownId = UUID.randomUUID().toString

  "CromwellApiService" should "return 404 for get of unknown workflow" in {
    Get(s"/workflows/${MockWorkflowManagerActor.unknownId}") ~>
      sealRoute(queryRoute) ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for get of a malformed workflow id" in {
    Get(s"/workflows/foobar") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return 200 for get of a known workflow id" in {
    Get(s"/workflows/${MockWorkflowManagerActor.runningWorkflowId}") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }

        assertResult(s"""{"id":"${MockWorkflowManagerActor.runningWorkflowId.toString}","status":"WorkflowRunning"}""") {
          responseAs[String]
        }
      }
  }
}
