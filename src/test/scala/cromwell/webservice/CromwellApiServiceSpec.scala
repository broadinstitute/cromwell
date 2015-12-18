package cromwell.webservice

import java.util.UUID

import akka.actor.{Actor, Props}
import akka.pattern.pipe
import cromwell.binding._
import cromwell.binding.values.{WdlFile, WdlInteger, WdlValue}
import cromwell.engine._
import cromwell.engine.backend.{CallLogs, WorkflowQueryResult}
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.MockWorkflowManagerActor.{submittedWorkflowId, unknownId}
import org.joda.time.DateTime
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import org.yaml.snakeyaml.Yaml
import spray.http.{DateTime => _, _}
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
  val submittedScatterWorkflowId = WorkflowId(UUID.randomUUID())
  val abortedWorkflowId = WorkflowId(UUID.randomUUID())

  def props: Props = Props(classOf[MockWorkflowManagerActor])
}

class MockWorkflowManagerActor extends Actor  {

  implicit lazy val hasher: WdlValue => SymbolHash = { x => SymbolHash("NOT REALLY IMPORTANT IN THESE TESTs!!") }

  val int8 = WdlInteger(8)
  val int8hash = int8.getHash

  val file = WdlFile("/tmp/ps.stdout.tmp")
  val fileHash = file.getHash

  def receive = {
    case SubmitWorkflow(sources) =>
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
            "three_step.cgrep.count" -> CallOutput(int8, int8hash),
            "three_step.ps.procs" -> CallOutput(file, fileHash),
            "three_step.wc.count" -> CallOutput(int8, int8hash)))
        case w => Future.failed(new WorkflowNotFoundException(s"Workflow '$w' not found"))
      }
      futureOutputs pipeTo sender

    case CallOutputs(id, callFqn) =>
      val futureOutputs =
        Future {
          id match {
            case MockWorkflowManagerActor.submittedWorkflowId =>
              callFqn match {
                case "three_step.cgrep" => Map("count" -> CallOutput(int8, int8hash))
                case "three_step.ps" => Map("procs" -> CallOutput(file, fileHash))
                case "three_step.wc" => Map("count" -> CallOutput(int8, int8hash))
                case _ => throw new CallNotFoundException(s"Bad call FQN: $callFqn")
              }
            case _ => throw new WorkflowNotFoundException(s"Bad workflow ID: $id")
          }
        }
      futureOutputs pipeTo sender

    case CallStdoutStderr(id, callFqn) =>
      val futureOutputs =
      Future {
        id match {
          case MockWorkflowManagerActor.submittedWorkflowId =>
            callFqn match {
              case "three_step.cgrep" => Seq(CallLogs(WdlFile("/path/to/cgrep-stdout"), WdlFile("/path/to/cgrep-stderr")))
              case "three_step.ps" => Seq(CallLogs(WdlFile("/path/to/ps-stdout"), WdlFile("/path/to/ps-stderr")))
              case "three_step.wc" => Seq(CallLogs(WdlFile("/path/to/wc-stdout"), WdlFile("/path/to/wc-stderr")))
              case "scatterwf.inside-scatter" =>
                Seq(CallLogs(WdlFile("/path/to/inside-scatter/shard0-stdout"), WdlFile("/path/to/inside-scatter/shard0-stderr")),
                    CallLogs(WdlFile("/path/to/inside-scatter/shard1-stdout"), WdlFile("/path/to/inside-scatter/shard1-stderr")))
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
            Map("three_step.ps" -> Seq(CallLogs(WdlFile("/path/to/ps-stdout"), WdlFile("/path/to/ps-stderr"))))
          case MockWorkflowManagerActor.submittedScatterWorkflowId =>
            Map("scatterwf.inside-scatter" ->
              Seq(CallLogs(WdlFile("/path/to/inside-scatter/shard0-stdout"), WdlFile("/path/to/inside-scatter/shard0-stderr")),
                  CallLogs(WdlFile("/path/to/inside-scatter/shard1-stdout"), WdlFile("/path/to/inside-scatter/shard1-stderr"))))
          case _ => throw new WorkflowNotFoundException(s"Bad workflow ID: $id")
        }
      }
      futureOutputs pipeTo sender

    case WorkflowQuery(rawParameters) =>
      val futureResult = Future {
        val head = rawParameters.head
        head match {
          case ("BadKey", _) =>
            // The exception text is rendered as the body, so there must be exception text or Spray will 500 (!)
            throw new IllegalArgumentException("aw snap")
          case ("status", _) =>
            WorkflowQueryResponse(
              Seq(
                WorkflowQueryResult(
                  id = UUID.randomUUID().toString,
                  name = "w",
                  status = "Succeeded",
                  start = new DateTime("2015-11-01T12:12:11"),
                  end = Option(new DateTime("2015-11-01T12:12:12")))
              ))
        }
      }
      futureResult pipeTo sender

    case CallCaching(id, parameters, callFqn) =>
      val parametersByKey = parameters.groupBy(_.key.toLowerCase.capitalize) mapValues { _ map { _.value } } mapValues { _.toSet }
      val futureResponse =
        Future {
          if (id == unknownId)
            throw new IllegalArgumentException("Unknown workflow")
          if (!parametersByKey.contains("Allow"))
            // Currently this is not strictly true as the "allow" validation only fails if "allow"s are non-boolean
            // or both true and false.  But really it would be better if "allow" was only specified once.
            throw new IllegalArgumentException("must specify 'allow' exactly once")
          if (parametersByKey.keys.size > 1)
            throw new IllegalArgumentException("Unrecognized parameters: " + (parametersByKey.keys.toSet - "allow").mkString(", "))
          if (callFqn.contains("bogus"))
            throw new IllegalArgumentException("Invalid call")
          // If we run the gauntlet of exception throwing checks, return a made up update count.
          1
        }
      futureResponse pipeTo sender
  }
}

class SwaggerServiceSpec extends FlatSpec with SwaggerService with ScalatestRouteTest with Matchers
with TableDrivenPropertyChecks {
  def actorRefFactory = system

  "Cromwell swagger docs" should "return 200" in {
    Get("/swagger/cromwell.yaml") ~>
      swaggerUiResourceRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult("2.0") {
          new Yaml()
            .loadAs(responseAs[String], classOf[java.util.Map[String, AnyRef]])
            .get("swagger")
        }
      }

    Get("/swagger/index.html") ~>
      swaggerUiResourceRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult("<!DOCTYPE html>") {
          responseAs[String].take(15)
        }
      }
  }

  "Cromwell swagger route" should "return 200" in {
    val pathExamples = Table("path", "/", "/swagger", "/swagger/cromwell.yaml", "/swagger/index.html", "/api",
      "/api/workflows/", "/api/workflows/v1", "/workflows/v1/outputs", "/workflows/v1/status",
      "/api/workflows/v1/validate", "/workflows", "/workflows/v1", "/workflows/v1/outputs", "/workflows/v1/status",
      "/workflows/v1/validate")

    forAll(pathExamples) { path =>
      Options(path) ~>
        swaggerUiResourceRoute ~>
        check {
          assertResult(StatusCodes.OK) {
            status
          }
          assertResult("OK") {
            responseAs[String]
          }
        }
    }
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
      sealRoute(statusRoute) ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for get of a malformed workflow id" in {
    Get(s"/workflows/$version/foobar/status") ~>
      statusRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return 200 for get of a known workflow id" in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.runningWorkflowId}/status") ~>
      statusRoute ~>
      check {
        assertResult(
          s"""{
              |  "id": "${MockWorkflowManagerActor.runningWorkflowId.toString}",
              |  "status": "Running"
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
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

  s"Cromwell submit workflow API $version" should "return 201 for a successful workflow submission " in {
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

  it should "return 400 for a malformed workflow inputs JSON " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> CromwellApiServiceSpec.MalformedInputsJson))) ~>
      submitRoute ~>
      check {
        assertResult("Expecting JSON object for workflowInputs and workflowOptions fields") {
          responseAs[String]
        }
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return 400 for a malformed workflow options JSON " in {
    Post("/workflows/$version", FormData(Seq("wdlSource" -> HelloWorld.wdlSource(), "workflowInputs" -> HelloWorld.rawInputs.toJson.toString(), "workflowOptions" -> CromwellApiServiceSpec.MalformedInputsJson))) ~>
      submitRoute ~>
      check {
        assertResult("Expecting JSON object for workflowInputs and workflowOptions fields") {
          responseAs[String]
        }
        assertResult(StatusCodes.BadRequest) {
          status
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
        assertResult(
          s"""{
              |  "valid": true
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.OK) {
          status
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
        assertResult(
          s"""{
              |  "valid": false,
              |  "error": "The following errors occurred while processing your inputs:\\n\\nRequired workflow input 'hello.hello.addressee' not specified."
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.BadRequest) {
          status
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
        assertResult(
          s"""{
              |  "valid": false,
              |  "error": "ERROR: Finished parsing without consuming all tokens.\\n\\nfoobar bad wdl!\\n^\\n     "
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.BadRequest) {
          status
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
        assertResult(
          s"""{
              |  "valid": false,
              |  "error": "Unexpected character 'o' at input index 0 (line 1, position 1), expected JSON Value:\\nfoobar bad json!\\n^\\n"
              |}""".stripMargin) {
          responseAs[String]
        }
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  s"Cromwell workflow outputs API $version" should "return 200 with GET of outputs on successful execution of workflow" in {
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

  "Cromwell call stdout/stderr API" should "return 200 with paths to stdout/stderr for calls inside a scatter" in {
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
  "Cromwell workflow stdout/stderr API" should "return 200 with paths to stdout/stderr with scattered calls" in {
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

  "Cromwell timings API" should "return 200 with an HTML document for the timings route"in {
    Get(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/timing") ~>
      timingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
        assertResult("<html>") {
          responseAs[String].substring(0, 6)
        }
      }
  }

  "Cromwell query API" should "return 400 for a bad query" in {
    Get(s"/workflows/$version/query?BadKey=foo") ~>
      queryRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }

  it should "return good results for a good query" in {
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

  "Cromwell single call caching API" should "work with good input" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call?allow=false") ~>
    callCachingRoute ~>
    check {
      assertResult(StatusCodes.OK) { status }
    }
  }

  it should "reject missing 'allow'" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call") ~>
    callCachingRoute ~>
    check {
      assertResult(StatusCodes.BadRequest) { status }
      assertResult(true) {
        responseAs[String].contains("must specify 'allow' exactly once")
      }
    }
  }

  it should "reject bogus calls" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/bogus?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Invalid call")
        }
      }
  }

  it should "reject invalid parameter keys" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching/w.good_call?allow=true&bogusKey=foo") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unrecognized parameters: ")
        }
      }
  }

  it should "reject bogus workflows" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.unknownId}/call-caching/w.good_call?allow=true") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unknown workflow")
        }
      }
  }

  "Cromwell all call caching API" should "work with good input" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching?allow=false") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.OK) { status }
      }
  }

  it should "reject missing 'allow'" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("must specify 'allow' exactly once")
        }
      }
  }

  it should "reject invalid parameter keys" in {
    Post(s"/workflows/$version/${MockWorkflowManagerActor.submittedScatterWorkflowId}/call-caching?allow=true&bogusKey=foo") ~>
      callCachingRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) { status }
        assertResult(true) {
          responseAs[String].contains("Unrecognized parameters: ")
        }
      }
  }

  it should "reject bogus workflows" in {
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
