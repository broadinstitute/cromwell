package cromwell.webservice

import java.util.UUID

import akka.actor.{Actor, Props}
import cromwell.core.{CallOutput, WorkflowId}
import cromwell.engine.backend.{CallLogs, WorkflowDescriptorBuilder, WorkflowQueryResult}
import cromwell.engine.workflow.OldStyleWorkflowManagerActor._
import cromwell.engine.{WorkflowAborted, WorkflowRunning}
import cromwell.server.WorkflowManagerSystem
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.MockWorkflowManagerActor.{submittedWorkflowId, unknownId}
import org.joda.time.DateTime
import org.scalatest.{FlatSpec, Matchers}
import spray.http.{DateTime => _, _}
import spray.json.DefaultJsonProtocol._
import spray.json._
import spray.testkit.ScalatestRouteTest
import wdl4s._
import wdl4s.values.{SymbolHash, WdlFile, WdlInteger, WdlValue}

import scala.util.{Failure, Success, Try}

object MockWorkflowManagerActor {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowRawInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerMessage

  val createdWorkflowId = WorkflowId(UUID.randomUUID())
  val runningWorkflowId = WorkflowId(UUID.randomUUID())
  val unknownId = WorkflowId(UUID.randomUUID())
  val submittedWorkflowId = WorkflowId(UUID.randomUUID())
  val submittedScatterWorkflowId = WorkflowId(UUID.randomUUID())
  val abortedWorkflowId = WorkflowId(UUID.randomUUID())

  def props: Props = Props(classOf[MockWorkflowManagerActor])
}

class MockWorkflowManagerActor extends Actor {
  this: WorkflowDescriptorBuilder =>

  implicit lazy val hasher: WdlValue => SymbolHash = { x => SymbolHash("NOT REALLY IMPORTANT IN THESE TESTs!!") }

  val int8 = WdlInteger(8)
  val int8hash = int8.computeHash

  val file = WdlFile("/tmp/ps.stdout.tmp")
  val fileHash = file.computeHash

  def receive = {
    case SubmitWorkflow(sources) =>
      val id = MockWorkflowManagerActor.submittedWorkflowId
      /*
        id will always succeed, but WorkflowDescriptor might not. If all of this works we
        want to pass Future[WorkflowId] but otherwise we want to pass the error. Try the
        WorkflowDescriptor - if it succeeds hand the id back in a future, otherwise the error
        from WorkflowDescriptor's validation
       */
      val message = Try(materializeWorkflowDescriptorFromSources(id = id, workflowSources = sources)) match {
        case Success(w) => WorkflowManagerSubmitSuccess(w.id)
        case Failure(e) => WorkflowManagerSubmitFailure(e)
      }
      sender ! message
    case WorkflowStatus(id) =>
      val message = id match {
        case MockWorkflowManagerActor.runningWorkflowId => WorkflowManagerStatusSuccess(id, WorkflowRunning)
        case MockWorkflowManagerActor.abortedWorkflowId => WorkflowManagerStatusSuccess(id, WorkflowAborted)
        case _ => WorkflowManagerStatusFailure(id, new WorkflowNotFoundException("Cromwell knows not this workflow"))
      }
      sender ! message
    case WorkflowAbort(id) =>
      val message = id match {
        case MockWorkflowManagerActor.runningWorkflowId => WorkflowManagerAbortSuccess(id)
        case MockWorkflowManagerActor.abortedWorkflowId =>
          WorkflowManagerAbortFailure(id, new IllegalStateException(s"Workflow ID '$id' is in terminal state 'Aborted' and cannot be aborted."))
        case x => WorkflowManagerAbortFailure(id, new WorkflowNotFoundException(s"Workflow '$x' not found."))
      }
      sender ! message
    case WorkflowOutputs(id) =>
      val message = id match {
        case MockWorkflowManagerActor.submittedWorkflowId =>
          WorkflowManagerWorkflowOutputsSuccess(id, Map(
            "three_step.cgrep.count" -> CallOutput(int8, Option(int8hash)),
            "three_step.ps.procs" -> CallOutput(file, Option(fileHash)),
            "three_step.wc.count" -> CallOutput(int8, Option(int8hash))))
        case w => WorkflowManagerWorkflowOutputsFailure(id, new WorkflowNotFoundException(s"Workflow '$w' not found"))
      }
      sender ! message
    case CallOutputs(id, callFqn) =>
      val message = id match {
        case MockWorkflowManagerActor.submittedWorkflowId =>
          callFqn match {
            case "three_step.cgrep" => WorkflowManagerCallOutputsSuccess(id, callFqn, Map("count" -> CallOutput(int8, Option(int8hash))))
            case "three_step.ps" => WorkflowManagerCallOutputsSuccess(id, callFqn, Map("procs" -> CallOutput(file, Option(fileHash))))
            case "three_step.wc" => WorkflowManagerCallOutputsSuccess(id, callFqn, Map("count" -> CallOutput(int8, Option(int8hash))))
            case _ => WorkflowManagerCallOutputsFailure(id, callFqn, new CallNotFoundException(s"Bad call FQN: $callFqn"))
          }
        case _ => WorkflowManagerCallOutputsFailure(id, callFqn, new WorkflowNotFoundException(s"Bad workflow ID: $id"))
      }
      sender ! message
    case CallStdoutStderr(id, callFqn) =>
      val message = id match {
        case MockWorkflowManagerActor.submittedWorkflowId =>
          callFqn match {
            case "three_step.cgrep" =>
              WorkflowManagerCallStdoutStderrSuccess(id, callFqn, Seq(CallLogs(WdlFile("/path/to/cgrep-stdout"), WdlFile("/path/to/cgrep-stderr"))))
            case "three_step.ps" =>
              WorkflowManagerCallStdoutStderrSuccess(id, callFqn, Seq(CallLogs(WdlFile("/path/to/ps-stdout"), WdlFile("/path/to/ps-stderr"))))
            case "three_step.wc" =>
              WorkflowManagerCallStdoutStderrSuccess(id, callFqn, Seq(CallLogs(WdlFile("/path/to/wc-stdout"), WdlFile("/path/to/wc-stderr"))))
            case "scatterwf.inside-scatter" =>
              WorkflowManagerCallStdoutStderrSuccess(id, callFqn, Seq(
                CallLogs(WdlFile("/path/to/inside-scatter/shard0-stdout"), WdlFile("/path/to/inside-scatter/shard0-stderr")),
                CallLogs(WdlFile("/path/to/inside-scatter/shard1-stdout"), WdlFile("/path/to/inside-scatter/shard1-stderr"))))
            case _ => WorkflowManagerCallStdoutStderrFailure(id, callFqn, new CallNotFoundException(s"Bad call FQN: $callFqn"))
          }
        case _ => WorkflowManagerCallStdoutStderrFailure(id, callFqn, new WorkflowNotFoundException(s"Bad workflow ID: $id"))
      }
      sender ! message
    case WorkflowStdoutStderr(id) =>
      val message = id match {
        case MockWorkflowManagerActor.submittedWorkflowId =>
          WorkflowManagerWorkflowStdoutStderrSuccess(id, Map("three_step.ps" -> Seq(CallLogs(WdlFile("/path/to/ps-stdout"), WdlFile("/path/to/ps-stderr")))))
        case MockWorkflowManagerActor.submittedScatterWorkflowId =>
          WorkflowManagerWorkflowStdoutStderrSuccess(id, Map("scatterwf.inside-scatter" ->
            Seq(CallLogs(WdlFile("/path/to/inside-scatter/shard0-stdout"), WdlFile("/path/to/inside-scatter/shard0-stderr")),
              CallLogs(WdlFile("/path/to/inside-scatter/shard1-stdout"), WdlFile("/path/to/inside-scatter/shard1-stderr")))))
        case _ => WorkflowManagerWorkflowStdoutStderrFailure(id, new WorkflowNotFoundException(s"Bad workflow ID: $id"))
      }
      sender ! message
    case WorkflowQuery(uri, rawParameters) =>
      val head = rawParameters.head
      head match {
        case ("BadKey", _) =>
          // The exception text is rendered as the body, so there must be exception text or Spray will 500 (!)
          sender ! WorkflowManagerQueryFailure(new IllegalArgumentException("Unrecognized query keys: BadKey"))
        case ("status", _) =>
          sender ! WorkflowManagerQuerySuccess(uri, WorkflowQueryResponse(
            Seq(
              WorkflowQueryResult(
                id = UUID.randomUUID().toString,
                name = "w",
                status = "Succeeded",
                start = new DateTime("2015-11-01T12:12:11"),
                end = Option(new DateTime("2015-11-01T12:12:12"))))),
            rawParameters.collectFirst { case (p, _) if p.contains("page") => QueryMetadata(Option(1), Option(5), Option(1)) })

      }
    case CallCaching(id, parameters, callFqn) =>
      val parametersByKey = parameters.groupBy(_.key.toLowerCase.capitalize) mapValues { _ map { _.value } } mapValues { _.toSet }
      val message =
        if (id == unknownId)
          WorkflowManagerCallCachingFailure(id, new IllegalArgumentException("Unknown workflow"))
        else if (!parametersByKey.contains("Allow"))
        // Currently this is not strictly true as the "allow" validation only fails if "allow"s are non-boolean
        // or both true and false.  But really it would be better if "allow" was only specified once.
          WorkflowManagerCallCachingFailure(id, new IllegalArgumentException("must specify 'allow' exactly once"))
        else if (parametersByKey.keys.size > 1)
          WorkflowManagerCallCachingFailure(id, new IllegalArgumentException("Unrecognized parameters: " + (parametersByKey.keys.toSet - "allow").mkString(", ")))
        else if (callFqn.contains("bogus"))
          WorkflowManagerCallCachingFailure(id, new IllegalArgumentException("Invalid call"))
        // If we run the gauntlet of checks, return a made up update count.
        else WorkflowManagerCallCachingSuccess(id, 1)
      sender ! message
  }
}

object CromwellApiServiceSpec {
  val MissingInputsJson: String = "{}"
  val MalformedInputsJson : String = "foobar bad json!"
  val MalformedWdl : String = "foobar bad wdl!"
}

class CromwellApiServiceSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {
  import spray.httpx.SprayJsonSupport._

  // BUG: Must be called once to statically initialize the backends, otherwise this Spec won't run if run alone.
  new WorkflowManagerSystem {}

  import akka.testkit._

  import scala.concurrent.duration._

  // The submit route takes a bit longer than the default 1 s while things initialize, when this spec is run by itself
  implicit val defaultTimeout = RouteTestTimeout(5.seconds.dilated)

  override def actorRefFactory = system
  override val workflowManager = actorRefFactory.actorOf(Props(new MockWorkflowManagerActor() with WorkflowDescriptorBuilder {
    override implicit  val actorSystem = context.system
  }))
  val version = "v1"

  // s"CromwellApiService $version" should "return 404 for get of unknown workflow" in {
  ignore should "return 404 for get of unknown workflow" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}") ~>
      sealRoute(statusRoute) ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  ignore should "return 400 for get of a malformed workflow id's status" in {
    Get(s"/workflows/$version/foobar/status") ~>
      statusRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          """{
             |  "status": "fail",
             |  "message": "Invalid workflow ID: 'foobar'."
             |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 200 for get of a known workflow id" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.runningWorkflowId}/status") ~>
      statusRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowManagerActor.runningWorkflowId.toString}",
              |  "status": "Running"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 404 for abort of unknown workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Workflow '${MockWorkflowManagerActor.unknownId}' not found."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 400 for abort of a malformed workflow id" in {
    Post(s"/workflows/$version/foobar/abort") ~>
      abortRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          """{
              |  "status": "fail",
              |  "message": "Invalid workflow ID: 'foobar'."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 403 for abort of a workflow in a terminal state" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.abortedWorkflowId}/abort") ~>
    abortRoute ~>
    check {
      assertResult(StatusCodes.Forbidden) {
        status
      }
      assertResult(
        s"""{
            |  "status": "error",
            |  "message": "Workflow ID '${MockWorkflowManagerActor.abortedWorkflowId}' is in terminal state 'Aborted' and cannot be aborted."
            |}""".stripMargin
      ) {
        responseAs[String]
      }
    }
  }

  ignore should "return 200 for abort of a known workflow id" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.runningWorkflowId}/abort") ~>
      abortRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowManagerActor.runningWorkflowId.toString}",
              |  "status": "Aborted"
              |}"""
            .stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 201 for a successful workflow submission " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> HelloWorld.rawInputs.toJson.toString()))) ~>
      submitRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowManagerActor.submittedWorkflowId.toString}",
              |  "status": "Submitted"
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.Created) {
          status
        }
      }
  }

  ignore should "return 400 for a malformed workflow inputs JSON " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> CromwellApiServiceSpec.MalformedInputsJson))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(true) {
          val fields: Map[String, JsValue] = responseAs[Map[String, JsValue]]
          fields.get("status").isDefined &&
          fields.get("status").get.asInstanceOf[JsString].value.equals("fail") &&
          fields.get("message").isDefined &&
          fields.get("message").get.asInstanceOf[JsString].value.contains("Workflow input processing failed.") &&
          fields.get("errors").isDefined &&
            (fields.get("errors").get match {
              case array: JsArray if array.elements.length == 1 =>
                array.elements.head.asInstanceOf[JsString].value.contains("contains invalid inputs JSON")
              case _ => false
            })
        }
      }
  }

  ignore should "return 400 for a malformed workflow options JSON " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> HelloWorld.rawInputs.toJson.toString(), "workflowOptions" -> CromwellApiServiceSpec.MalformedInputsJson))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        val fields: Map[String, JsValue] = responseAs[Map[String, JsValue]]
        assertResult(true) {
          fields.get("status").isDefined &&
            fields.get("status").get.asInstanceOf[JsString].value.equals("fail")
        }
        assertResult(true) {
          fields.get("message").isDefined &&
            fields.get("message").get.asInstanceOf[JsString].value.contains("Workflow input processing failed")
        }
        assertResult(true) {
          fields.get("errors").isDefined &&
            (fields.get("errors").get match {
              case array: JsArray if array.elements.length == 1 =>
                array.elements.head.asInstanceOf[JsString].value.contains("contains invalid options JSON")
              case _ => false
            })
        }
      }
  }

  ignore should "return 200 for a successful workflow submission " in {
    val inputs = HelloWorld.rawInputs.toJson

    Post("/workflows/$version/batch",
      FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> s"[$inputs, $inputs]"))) ~>
      submitBatchRoute ~>
      check {
        assertResult(
          s"""[{
              |  "id": "${MockWorkflowManagerActor.submittedWorkflowId.toString}",
              |  "status": "Submitted"
              |}, {
              |  "id": "${MockWorkflowManagerActor.submittedWorkflowId.toString}",
              |  "status": "Submitted"
              |}]""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  // TODO: Test tha batch submission returns expected workflow ids in order
  // TODO: Also (assuming we still validate on submit) test a batch of mixed inputs that return submitted and failed

  ignore should "return 200 with GET of outputs on successful execution of workflow" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.submittedWorkflowId.toString}/outputs") ~>
      workflowOutputsRoute ~>
      check {
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
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 404 with outputs on unknown workflow" in {
    Get(s"/workflows/$version/$unknownId/outputs") ~>
    workflowOutputsRoute ~>
    check {
      assertResult(StatusCodes.NotFound) {
        status
      }
      assertResult(
        s"""{
            |  "status": "error",
            |  "message": "Workflow '$unknownId' not found."
            |}""".stripMargin
      ) {
        responseAs[String]
      }
    }
  }

  ignore should "return 200 with outputs on successful execution of workflow" in {
    Get(s"/workflows/$version/$submittedWorkflowId/outputs/three_step.wc") ~>
    callOutputsRoute ~>
    check {
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
      assertResult(StatusCodes.OK) {
        status
      }
    }
  }

  ignore should "return 404 for unknown workflow" in {
    Get(s"/workflows/$version/$unknownId/outputs/three_step.wc") ~>
      callOutputsRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Workflow '$unknownId' not found."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 404 for unknown call ID" in {
    Get(s"/workflows/$version/$submittedWorkflowId/outputs/bogus_workflow.bogus_call") ~>
      callOutputsRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Call bogus_workflow.bogus_call not found for workflow '$submittedWorkflowId'."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 400 for malformed workflow ID" in {
    Get(s"/workflows/$version/foobar/outputs/three_step.wc") ~>
      callOutputsRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Invalid workflow ID: 'foobar'."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 405 with POST of outputs on successful execution of workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedWorkflowId.toString}/outputs") ~>
      sealRoute(workflowOutputsRoute) ~>
      check {
        assertResult(StatusCodes.MethodNotAllowed) {
          status
        }
      }
  }

  ignore should "return 200 with paths to stdout/stderr for a call" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs/three_step.wc") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "$submittedWorkflowId",
              |  "logs": {
              |    "three_step.wc": [{
              |      "stdout": "/path/to/wc-stdout",
              |      "stderr": "/path/to/wc-stderr"
              |    }]
              |  }
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 200 with paths to stdout/stderr for calls inside a scatter" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs/scatterwf.inside-scatter") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "$submittedWorkflowId",
              |  "logs": {
              |    "scatterwf.inside-scatter": [{
              |      "stdout": "/path/to/inside-scatter/shard0-stdout",
              |      "stderr": "/path/to/inside-scatter/shard0-stderr"
              |    }, {
              |      "stdout": "/path/to/inside-scatter/shard1-stdout",
              |      "stderr": "/path/to/inside-scatter/shard1-stderr"
              |    }]
              |  }
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 404 if the workflow ID is not found" in {
    val randomID = {UUID.randomUUID().toString}
    Get(s"/workflows/$version/$randomID/logs/three_step.wc") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Workflow '$randomID' not found."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 404 if the Call FQN is not found" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs/three_step.wcBADBAD") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
        assertResult(
          s"""{
              |  "status": "error",
              |  "message": "Call three_step.wcBADBAD not found for workflow '$submittedWorkflowId'."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 400 for get of a malformed workflow id's logs" in {
    Get(s"/workflows/$version/foobar/logs/three_step.wc") ~>
      callStdoutStderrRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Invalid workflow ID: 'foobar'."
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return 200 with paths to stdout/stderr" in {
    Get(s"/workflows/$version/$submittedWorkflowId/logs") ~>
      workflowStdoutStderrRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "$submittedWorkflowId",
              |  "logs": {
              |    "three_step.ps": [{
              |      "stdout": "/path/to/ps-stdout",
              |      "stderr": "/path/to/ps-stderr"
              |    }]
              |  }
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 200 with paths to stdout/stderr with scattered calls" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/logs") ~>
      workflowStdoutStderrRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowManagerActor.submittedScatterWorkflowId}",
              |  "logs": {
              |    "scatterwf.inside-scatter": [{
              |      "stdout": "/path/to/inside-scatter/shard0-stdout",
              |      "stderr": "/path/to/inside-scatter/shard0-stderr"
              |    }, {
              |      "stdout": "/path/to/inside-scatter/shard1-stdout",
              |      "stderr": "/path/to/inside-scatter/shard1-stderr"
              |    }]
              |  }
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  ignore should "return 200 with an HTML document for the timings route"in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/timing") ~>
      timingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
        assertResult("<html>") {
          responseAs[String].substring(0, 6)
        }
      }
  }

  ignore should "return 400 for a bad query" in {
    Get(s"/workflows/$version/query?BadKey=foo") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Unrecognized query keys: BadKey"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return good results for a good query" in {
    Get(s"/workflows/$version/query?status=Succeeded") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\",")
        }
      }
  }

  it should "return link headers for pagination when page and pagesize are set for a good query" in {
    Get(s"/workflows/$version/query?status=Succeeded&page=1&pagesize=5") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\",")
          (headers count { header => header.is("link") }) == 4
        }
      }
  }

  behavior of "Cromwell query post API"

  ignore should "return 400 for a bad query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"BadKey":"foo"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult(
          s"""{
              |  "status": "fail",
              |  "message": "Unrecognized query keys: BadKey"
              |}""".stripMargin
        ) {
          responseAs[String]
        }
      }
  }

  ignore should "return good results for a good query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\",")
        }
      }
  }

  it should "return link headers for pagination when page and pagesize are set for a good query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":"Succeeded"},  {"page": "1"}, {"pagesize": "5"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\",")
          (headers count { header => header.is("link") }) == 4
        }
      }
  }


  ignore should "return good results for a multiple query map body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`,
      """[{"status":"Succeeded"}, {"status":"Failed"}]""")) ~>
      queryPostRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult(true) {
          body.asString.contains("\"status\": \"Succeeded\",")
        }
      }
  }

  ignore should "return 400 bad request for a bad query format body" in {
    Post(s"/workflows/$version/query", HttpEntity(ContentTypes.`application/json`, """[{"status":["Succeeded"]}]""")) ~>
      sealRoute(queryPostRoute) ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  behavior of "Cromwell single call caching API"

  ignore should "disallow call caching for a call" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call?allow=false") ~>
    callCachingRoute ~>
    check {
      assertResult(StatusCodes.OK) { status }
    }
  }

  ignore should "reject missing 'allow' when disabling call caching for a call" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call") ~>
    callCachingRoute ~>
    check {
      assertResult(StatusCodes.BadRequest) { status }
      assertResult(true) {
        responseAs[String].contains("must specify 'allow' exactly once")
      }
    }
  }

  ignore should "reject bogus calls" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/bogus?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Invalid call")
        }
      }
  }

  ignore should "reject invalid parameter keys when enabling call caching for a call" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call?allow=true&bogusKey=foo") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unrecognized parameters: ")
        }
      }
  }

  ignore should "reject bogus workflows when enabling call caching for a call" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/call-caching/w.good_call?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unknown workflow")
        }
      }
  }

  ignore should "disallow call caching for a workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching?allow=false") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
      }
  }

  ignore should "reject missing 'allow' when disabling call caching for a workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("must specify 'allow' exactly once")
        }
      }
  }

  ignore should "reject invalid parameter keys when enabling call caching for a workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching?allow=true&bogusKey=foo") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unrecognized parameters: ")
        }
      }
  }

  ignore should "reject bogus workflows when enabling call caching for a workflow" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/call-caching?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unknown workflow")
        }
      }
  }
}
