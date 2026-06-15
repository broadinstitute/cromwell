package cromwell.webservice.routes.wes

import akka.actor.Props
import akka.http.scaladsl.testkit.ScalatestRouteTest
import akka.util.Timeout
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.webservice.routes.CromwellApiService
import cromwell.webservice.routes.CromwellApiServiceSpec.MockWorkflowStoreActor
import org.scalatest.flatspec.AsyncFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.concurrent.ExecutionContextExecutor

class ServiceInfoSpec extends AsyncFlatSpec with ScalatestRouteTest with Matchers {
  val actorRefFactory = system

  // ExecutionContextExecutor rather than ExecutionContext to make sure we take precedence
  // over the `executor: ExecutionContextExecutor` inherited from ScalatestRouteTest
  implicit val ec: ExecutionContextExecutor = system.dispatcher
  implicit val timeout: Timeout = 5.seconds

  val workflowStoreActor = actorRefFactory.actorOf(Props(new MockWorkflowStoreActor()))

  CromwellLanguages.initLanguages(LanguageConfiguration.AllLanguageEntries)

  behavior of "ServiceInfo"

  val expectedResponse = WesStatusInfoResponse(
    Map("WDL" -> Set("1.1", "1.0", "draft-2", "cascades")),
    List("1.0"),
    Set("ftp", "s3", "drs", "gcs", "http"),
    Map("Cromwell" -> CromwellApiService.cromwellVersion),
    List(),
    Map(WesState.Running -> 5, WesState.Queued -> 3, WesState.Canceling -> 2),
    "https://cromwell.readthedocs.io/en/stable/",
    "https://cromwell.readthedocs.io/en/stable/",
    Map()
  )

  it should "should eventually build the right WesResponse" in {
    ServiceInfo.toWesResponse(workflowStoreActor) map { r =>
      assert(r == expectedResponse)
    }
  }
}
