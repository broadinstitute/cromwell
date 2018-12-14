package cromwell.webservice.routes.wes

import akka.actor.Props
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.model.HttpMethods.POST
import akka.http.scaladsl.testkit.ScalatestRouteTest
import cromwell.webservice.routes.CromwellApiServiceSpec

import scala.concurrent.duration._
import cromwell.webservice.routes.CromwellApiServiceSpec.{MockServiceRegistryActor, MockWorkflowManagerActor, MockWorkflowStoreActor}
import org.scalatest.{AsyncFlatSpec, Matchers}
import WesResponseJsonSupport._
import akka.http.scaladsl.server.MethodRejection
import spray.json._

class WesRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers with WesRouteSupport {
  val actorRefFactory = system
  override implicit val ec = system.dispatcher
  override implicit val timeout = 5.seconds

  override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
  override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
  override val metadataBuilderRegulatorActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor())) // Not using a regulator for tests
  override val workflowManagerActor = actorRefFactory.actorOf(Props(new MockWorkflowManagerActor()))

  val version = "v1"

  behavior of "WES API /status endpoint"
  it should "return PAUSED when on hold" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.OnHoldWorkflowId}/status") ~>
      wesRoutes ~>
        check {
          responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.OnHoldWorkflowId.toString, WesState.Paused)
        }
  }

  it should "return QUEUED when submitted" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.ExistingWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.ExistingWorkflowId.toString, WesState.Queued)
      }
  }

  it should "return RUNNING when running" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.RunningWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.RunningWorkflowId.toString, WesState.Running)
      }
  }

  it should "return CANCELING when aborting" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.AbortingWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.AbortingWorkflowId.toString, WesState.Canceling)
      }
  }

  it should "return CANCELED when on hold" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.AbortedWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.AbortedWorkflowId.toString, WesState.Canceled)
      }
  }

  it should "return COMPLETE when successful" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.SucceededWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.SucceededWorkflowId.toString, WesState.Complete)
      }
  }

  it should "return EXECUTOR_ERROR when a workflow fails" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.FailedWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.FailedWorkflowId.toString, WesState.ExecutorError)
      }
  }

  behavior of "WES API /{run_id} endpoint"
  it should "return 404 when requesting an unknown workflow" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.UnrecognizedWorkflowId}") ~>
      wesRoutes ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }

        responseAs[WesErrorResponse] shouldEqual WesErrorResponse("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)
      }
  }

  it should "correctly return a WES run log for a known workflow" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.SupWdlWorkflowId}") ~>
      wesRoutes ~>
      check {
        responseAs[JsObject] shouldEqual WesRouteSupportSpec.Foo
      }
  }

  behavior of "WES API /cancel endpoint"
  it should "return 404 for abort of unknown workflow" in {
    val workflowId = CromwellApiServiceSpec.UnrecognizedWorkflowId

    Post(s"/ga4gh/wes/$version/runs/$workflowId/cancel") ~>
      wesRoutes ~>
      check {
        assertResult(StatusCodes.NotFound) {
          status
        }
      }
  }

  it should "return 400 for abort of a malformed workflow id" in {
    Post(s"/ga4gh/wes/$version/runs/foobar/cancel") ~>
      wesRoutes ~>
      check {
        assertResult(StatusCodes.InternalServerError) {
          status
        }

        responseAs[WesErrorResponse] shouldEqual WesErrorResponse("Invalid workflow ID: 'foobar'.", StatusCodes.InternalServerError.intValue)
      }
  }

  it should "return 200 for abort of an OnHold workflow" in {
    Post(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.OnHoldWorkflowId}/cancel") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunId] shouldEqual WesRunId(CromwellApiServiceSpec.OnHoldWorkflowId.toString)

        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  it should "return 200 for abort of a workflow in Submitted state" in {
    Post(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.SubmittedWorkflowId}/cancel") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunId] shouldEqual WesRunId(CromwellApiServiceSpec.SubmittedWorkflowId.toString)

        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  it should "return 200 for abort of an aborting workflow" in {
    Post(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.AbortingWorkflowId}/cancel") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunId] shouldEqual WesRunId(CromwellApiServiceSpec.AbortingWorkflowId.toString)

        assertResult(StatusCodes.OK) {
          status
        }
      }
  }

  it should "fail if you submit via GET" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.ExistingWorkflowId}/cancel") ~>
      wesRoutes ~>
      check {
        rejection shouldEqual MethodRejection(POST)
      }
  }
}

object WesRouteSupportSpec {
  val Foo = """{
              |  "request": {
              |    "workflow_url": "",
              |    "workflow_engine_parameters": {},
              |    "workflow_type": "WDL",
              |    "workflow_params": {
              |      "supsup.sup.addressee": "dog"
              |    },
              |    "tags": {},
              |    "workflow_type_version": "draft-2"
              |  },
              |  "state": "COMPLETE",
              |  "outputs": {
              |    "supsup.sup.salutation": "yo sup dog!"
              |  },
              |  "run_id": "cab8b246-8c78-4ca9-b2c3-1c8377da986e",
              |  "task_logs": [
              |    {
              |      "name": "supsup.sup",
              |      "start_time": "2018-12-14T14:25:09.376-05:00",
              |      "stdout": "/Users/jgentry/projects/cromwell/cromwell-executions/supsup/cab8b246-8c78-4ca9-b2c3-1c8377da986e/call-sup/execution/stdout",
              |      "end_time": "2018-12-14T14:25:24.530-05:00",
              |      "exit_code": 0,
              |      "cmd": [
              |        "sleep 10\necho \"yo sup dog!\""
              |      ],
              |      "stderr": "/Users/jgentry/projects/cromwell/cromwell-executions/supsup/cab8b246-8c78-4ca9-b2c3-1c8377da986e/call-sup/execution/stderr"
              |    }
              |  ],
              |  "run_log": {
              |    "name": "supsup",
              |    "start_time": "2018-12-14T14:25:08.124-05:00",
              |    "end_time": "2018-12-14T14:25:25.703-05:00"
              |  }
              |}""".stripMargin.parseJson.asInstanceOf[JsObject]
}