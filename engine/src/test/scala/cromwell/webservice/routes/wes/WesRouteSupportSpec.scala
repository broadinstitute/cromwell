package cromwell.webservice.routes.wes

import akka.actor.Props
import akka.http.scaladsl.model.HttpMethods.POST
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.MethodRejection
import akka.http.scaladsl.testkit.{RouteTestTimeout, ScalatestRouteTest}
import cromwell.util.SampleWdl.HelloWorld
import cromwell.webservice.routes.CromwellApiServiceSpec
import cromwell.webservice.routes.CromwellApiServiceSpec.{
  MockServiceRegistryActor,
  MockWorkflowManagerActor,
  MockWorkflowStoreActor
}
import cromwell.webservice.routes.wes.WesResponseJsonSupport._
import org.scalatest.flatspec.AsyncFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

import scala.concurrent.duration._

class WesRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers with WesRouteSupport {

  val actorRefFactory = system
  implicit override val ec = system.dispatcher
  override val timeout = routeTestTimeout.duration
  implicit def routeTestTimeout = RouteTestTimeout(5.seconds)

  override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
  override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
  override val workflowManagerActor = actorRefFactory.actorOf(Props(new MockWorkflowManagerActor()))

  val version = "v1"

  behavior of "WES API /status endpoint"
  it should "return PAUSED when on hold" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.OnHoldWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.OnHoldWorkflowId.toString,
                                                          WesState.Paused
        )
      }
  }

  it should "return QUEUED when submitted" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.ExistingWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.ExistingWorkflowId.toString,
                                                          WesState.Queued
        )
      }
  }

  it should "return RUNNING when running" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.RunningWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.RunningWorkflowId.toString,
                                                          WesState.Running
        )
      }
  }

  it should "return CANCELING when aborting" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.AbortingWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.AbortingWorkflowId.toString,
                                                          WesState.Canceling
        )
      }
  }

  it should "return CANCELED when on hold" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.AbortedWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.AbortedWorkflowId.toString,
                                                          WesState.Canceled
        )
      }
  }

  it should "return COMPLETE when successful" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.SucceededWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.SucceededWorkflowId.toString,
                                                          WesState.Complete
        )
      }
  }

  it should "return EXECUTOR_ERROR when a workflow fails" in {
    Get(s"/ga4gh/wes/$version/runs/${CromwellApiServiceSpec.FailedWorkflowId}/status") ~>
      wesRoutes ~>
      check {
        responseAs[WesRunStatus] shouldEqual WesRunStatus(CromwellApiServiceSpec.FailedWorkflowId.toString,
                                                          WesState.ExecutorError
        )
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

        responseAs[WesErrorResponse] shouldEqual WesErrorResponse("Invalid workflow ID: 'foobar'.",
                                                                  StatusCodes.InternalServerError.intValue
        )
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

  behavior of "WES API /runs POST endpoint"
  it should "return 201 for a successful workflow submission" in {
    val workflowSource = Multipart.FormData.BodyPart(
      "workflow_url",
      HttpEntity(
        MediaTypes.`application/json`,
        "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"
      )
    )
    val workflowInputs =
      Multipart.FormData.BodyPart("workflow_params",
                                  HttpEntity(MediaTypes.`application/json`, HelloWorld.rawInputs.toJson.toString())
      )
    val formData = Multipart.FormData(workflowSource, workflowInputs).toEntity()
    Post(s"/ga4gh/wes/$version/runs", formData) ~>
      wesRoutes ~>
      check {
        assertResult(s"""{
                        |  "run_id": "${CromwellApiServiceSpec.ExistingWorkflowId.toString}"
                        |}""".stripMargin) {
          responseAs[String].parseJson.prettyPrint
        }
        assertResult(StatusCodes.Created) {
          status
        }
        headers should be(Seq.empty)
      }
  }

  behavior of "WES API /runs GET endpoint"
  it should "return results for a good query" in {
    Get(s"/ga4gh/wes/v1/runs") ~>
      wesRoutes ~>
      check {
        status should be(StatusCodes.OK)
        contentType should be(ContentTypes.`application/json`)
        val results = responseAs[JsObject].fields("runs").convertTo[Seq[JsObject]]
        results.head.fields("run_id") should be(JsString(CromwellApiServiceSpec.ExistingWorkflowId.toString))
        results.head.fields("state") should be(JsString("COMPLETE"))
      }
  }

  behavior of "WES API /runs/{run_id} endpoint"
  it should "return valid metadata when supplied a run_id" in {
    Get(s"/ga4gh/wes/v1/runs/${CromwellApiServiceSpec.wesWorkflowId}") ~>
      wesRoutes ~>
      check {
        status should be(StatusCodes.OK)
        val result = responseAs[JsObject]
        result.fields.keys should contain allOf ("request", "run_id", "state")
        result.fields("state") should be(JsString("RUNNING"))
        result.fields("run_id") should be(JsString(CromwellApiServiceSpec.wesWorkflowId.toString))
      }
  }
}
