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

class WesRouteSupportSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers with WesRouteSupport {
  val actorRefFactory = system
  override implicit val ec = system.dispatcher
  override implicit val timeout = 5.seconds

  override val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))
  override val serviceRegistryActor = actorRefFactory.actorOf(Props(new MockServiceRegistryActor()))
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
