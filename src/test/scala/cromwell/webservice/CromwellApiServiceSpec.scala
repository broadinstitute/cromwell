package cromwell.webservice

import java.util.UUID

import akka.actor.{Actor, Props}
import akka.pattern.pipe
import cromwell.binding._
import cromwell.binding.values.{WdlFile, WdlInteger}
import cromwell.engine._
import cromwell.engine.backend.StdoutStderr
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.MockWorkflowManagerActor.{submittedWorkflowId, unknownId}
import org.scalatest.{FlatSpec, Matchers}
import spray.http._
import spray.json.DefaultJsonProtocol._
import spray.json._
import spray.testkit.ScalatestRouteTest

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future

object MockWorkflowManagerActor {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowRawInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerMessage

  val createdWorkflowId = WorkflowId(UUID.randomUUID())
  val runningWorkflowId = WorkflowId(UUID.randomUUID())
  val unknownId = WorkflowId(UUID.randomUUID())
  val submittedWorkflowId = WorkflowId(UUID.randomUUID())
  val abortedWorkflowId = WorkflowId(UUID.randomUUID())

  def props: Props = Props(classOf[MockWorkflowManagerActor])
}

class MockWorkflowManagerActor extends Actor  {

  def receive = {
    case SubmitWorkflow(wdlSource, wdlJson, rawInputs) =>
      sender ! MockWorkflowManagerActor.submittedWorkflowId

    case WorkflowStatus(id) =>
      val msg = id match {
        case MockWorkflowManagerActor.runningWorkflowId =>
          Some(WorkflowRunning)
        case MockWorkflowManagerActor.abortedWorkflowId =>
          Some(WorkflowAborted)
        case _ =>
          None
      }
      sender ! msg

    case WorkflowAbort(id) =>
      val msg = id match {
        case MockWorkflowManagerActor.runningWorkflowId =>
          Some(WorkflowRunning)
        case _ =>
          None
      }
      sender ! msg

    case WorkflowOutputs(id) =>
      val futureOutputs = id match {
        case MockWorkflowManagerActor.submittedWorkflowId =>
          Future.successful(Map(
            "three_step.cgrep.count" -> WdlInteger(8),
            "three_step.ps.procs" -> WdlFile("/tmp/ps.stdout.tmp"),
            "three_step.wc.count" -> WdlInteger(8)))
        case w => Future.failed(new WorkflowNotFoundException(s"Workflow '$w' not found"))
      }
      futureOutputs pipeTo sender

    case CallOutputs(id, callFqn, _) =>
      val futureOutputs =
        Future {
          id match {
            case MockWorkflowManagerActor.submittedWorkflowId =>
              callFqn match {
                case "three_step.cgrep" => Map("count" -> WdlInteger(8))
                case "three_step.ps" => Map("procs" -> WdlFile("/tmp/ps.stdout.tmp"))
                case "three_step.wc" => Map("count" -> WdlInteger(8))
                case _ => throw new CallNotFoundException(s"Bad call FQN: $callFqn")
              }
            case _ => throw new WorkflowNotFoundException(s"Bad workflow ID: $id")
          }
        }
      futureOutputs pipeTo sender

    case CallStdoutStderr(id, callFqn, _) =>
      val futureOutputs =
      Future {
        id match {
          case MockWorkflowManagerActor.submittedWorkflowId =>
            callFqn match {
              case "three_step.cgrep" => StdoutStderr(WdlFile("/path/to/cgrep-stdout"), WdlFile("/path/to/cgrep-stderr"))
              case "three_step.ps" => StdoutStderr(WdlFile("/path/to/ps-stdout"), WdlFile("/path/to/ps-stderr"))
              case "three_step.wc" => StdoutStderr(WdlFile("/path/to/wc-stdout"), WdlFile("/path/to/wc-stderr"))
              case _ => throw new CallNotFoundException(s"Bad call FQN: $callFqn")
            }
          case _ => throw new WorkflowNotFoundException(s"Bad workflow ID: $id")
        }
      }
      futureOutputs pipeTo sender

    case WorkflowStdoutStderr(id) =>
      val futureOutputs =
      Future {
        id match {
          case MockWorkflowManagerActor.submittedWorkflowId =>
            Map("three_step.ps" -> StdoutStderr(WdlFile("/path/to/ps-stdout"), WdlFile("/path/to/ps-stderr")))
          case _ => throw new WorkflowNotFoundException(s"Bad workflow ID: $id")
        }
      }
      futureOutputs pipeTo sender
  }
}

object CromwellApiServiceSpec {
  val MissingInputsJson: String = "{}"
  val MalformedInputsJson : String = "foobar bad json!"
  val MalformedWdl : String = "foobar bad wdl!"
}

class CromwellApiServiceSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {
  def actorRefFactory = system
  val workflowManager = system.actorOf(Props(new MockWorkflowManagerActor()))
  val version = "v1"

  s"CromwellApiService $version" should "return 404 for get of unknown workflow" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}") ~>
      sealRoute(queryRoute) ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for get of a malformed workflow id" in {
    Get(s"/workflows/$version/foobar/status") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return 200 for get of a known workflow id" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.runningWorkflowId}/status") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }

        assertResult(
          s"""{
             |  "id": "${MockWorkflowManagerActor.runningWorkflowId.toString}",
             |  "status": "Running"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  "CromwellApiService" should "return 404 for abort of unknown workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for abort of a malformed workflow id" in {
    Post(s"/workflows/$version/foobar/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return 403 for abort of a workflow in a terminal state" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.abortedWorkflowId}/abort") ~>
    abortRoute ~>
    check {
      assertResult(StatusCodes.Forbidden) {
        status
      }
    }
  }

  it should "return 200 for abort of a known workflow id" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.runningWorkflowId}/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }

        assertResult(
          s"""{
             |  "id": "${MockWorkflowManagerActor.runningWorkflowId.toString}",
             |  "status": "Aborted"
             |}"""
            .stripMargin) {
          responseAs[String]
        }
      }
  }

  s"Cromwell submit workflow API $version" should "return 201 for a successful workflow submission " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> HelloWorld.rawInputs.toJson.toString()))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.Created) {
          status
        }
        assertResult(
          s"""{
             |  "id": "${MockWorkflowManagerActor.submittedWorkflowId.toString}",
             |  "status": "Submitted"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return 400 for a malformed JSON " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> CromwellApiServiceSpec.MalformedInputsJson))) ~>
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

  s"Cromwell validate workflow API $version" should "return 200 for a successful workflow submission" in {
    Post("/workflows/$version/validate", FormData(Seq(
      "wdlSource" -> HelloWorld.wdlSource(),
      "workflowInputs" -> HelloWorld.rawInputs.toJson.toString()))
    ) ~>
      validateRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(
          s"""{
             |  "valid": true
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return 400 for a missing input" in {
    Post("/workflows/$version/validate", FormData(Seq(
      "wdlSource" -> HelloWorld.wdlSource(),
      "workflowInputs" -> CromwellApiServiceSpec.MissingInputsJson))
    ) ~>
      validateRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
             |  "valid": false,
             |  "error": "The following errors occurred while processing your inputs:\\n\\nRequired workflow input 'hello.hello.addressee' not specified."
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return 400 for a malformed WDL" in {
    Post(s"/workflows/$version/validate", FormData(Seq(
      "wdlSource" -> CromwellApiServiceSpec.MalformedWdl,
      "workflowInputs" -> HelloWorld.rawInputs.toJson.toString()))
    ) ~>
      validateRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
             |  "valid": false,
             |  "error": "ERROR: Finished parsing without consuming all tokens.\\n\\nfoobar bad wdl!\\n^\\n     "
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return 400 for a malformed JSON" in {
    Post("/workflows/$version/validate", FormData(Seq(
      "wdlSource" -> HelloWorld.wdlSource(),
      "workflowInputs" -> CromwellApiServiceSpec.MalformedInputsJson))
    ) ~>
      validateRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
             |  "valid": false,
             |  "error": "Unexpected character 'o' at input index 0 (line 1, position 1), expected JSON Value:\\nfoobar bad json!\\n^\\n"
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  s"Cromwell workflow outputs API $version" should "return 200 with GET of outputs on successful execution of workflow" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.submittedWorkflowId.toString}/outputs") ~>
      workflowOutputsRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(
          s"""{
             |  "id": "${MockWorkflowManagerActor.submittedWorkflowId.toString}",
             |  "outputs": {
             |    "three_step.cgrep.count": 8,
             |    "three_step.ps.procs": "/tmp/ps.stdout.tmp",
             |    "three_step.wc.count": 8
             |  }
             |}""".stripMargin) {
            responseAs[String]
          }
      }
  }

  it should "return 404 with outputs on unknown workflow" in {
    Get(s"/workflows/$version/$unknownId/outputs") ~>
    workflowOutputsRoute ~>
    check {
      assertResult(StatusCodes.NotFound) {
        status
      }
    }
  }

  "Cromwell call outputs API" should "return 200 with outputs on successful execution of workflow" in {
    Get(s"/workflows/$version/$submittedWorkflowId/outputs/three_step.wc") ~>
    callOutputsRoute ~>
    check {
      assertResult(StatusCodes.OK) {
        status
      }
      assertResult(
        s"""{
           |  "id": "$submittedWorkflowId",
           |  "callFqn": "three_step.wc",
           |  "outputs": {
           |    "count": 8
           |  }
           |}""".stripMargin) {
        responseAs[String]
      }
    }
  }

  it should "return 404 for unknown workflow" in {
    Get(s"/workflows/$version/$unknownId/outputs/three_step.wc") ~>
      callOutputsRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 404 for unknown call ID" in {
    Get(s"/workflows/$version/$submittedWorkflowId/outputs/bogus_workflow.bogus_call") ~>
      callOutputsRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for malformed workflow ID" in {
    Get(s"/workflows/$version/foobar/outputs/three_step.wc") ~>
      callOutputsRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }
  it should "return 405 with POST of outputs on successful execution of workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedWorkflowId.toString}/outputs") ~>
      sealRoute(workflowOutputsRoute) ~>
      check {
        assertResult(StatusCodes.MethodNotAllowed) {
          status
        }
      }
  }

  "Cromwell call stdout/stderr API" should "return 200 with paths to stdout/stderr" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs/three_step.wc") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(
          s"""{
             |  "id": "$submittedWorkflowId",
             |  "logs": {
             |    "three_step.wc": {
             |      "stdout": "/path/to/wc-stdout",
             |      "stderr": "/path/to/wc-stderr"
             |    }
             |  }
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

  it should "return 404 if the workflow ID is not found" in {
    Get(s"/workflows/$version/${UUID.randomUUID().toString}/logs/three_step.wc") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 404 if the Call FQN is not found" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs/three_step.wcBADBAD") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for get of a malformed workflow id" in {
    Get(s"/workflows/$version/foobar/logs/three_step.wc") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  "Cromwell workflow stdout/stderr API" should "return 200 with paths to stdout/stderr" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs") ~>
      workflowStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(
          s"""{
             |  "id": "$submittedWorkflowId",
             |  "logs": {
             |    "three_step.ps": {
             |      "stdout": "/path/to/ps-stdout",
             |      "stderr": "/path/to/ps-stderr"
             |    }
             |  }
             |}""".stripMargin) {
          responseAs[String]
        }
      }
  }

}
