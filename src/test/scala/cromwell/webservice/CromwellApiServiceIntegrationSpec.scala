package cromwell.webservice

import cromwell.CromwellSpec
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.{FlatSpec, Matchers}
import spray.http.{FormData, StatusCodes}
import spray.json.DefaultJsonProtocol._
import spray.json._
import spray.testkit.ScalatestRouteTest

import scala.concurrent.Await
import scala.concurrent.duration.Duration


class CromwellApiServiceIntegrationSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {
  def actorRefFactory = system
  val dataAccess = DataAccess.instance
  val workflowManager = system.actorOf(WorkflowManagerActor.props(CromwellSpec.BackendInstance))
  val version = "v1"

  override protected def afterAll() {
    super.afterAll()
    Await.result(dataAccess.shutdown(), Duration.Inf)
  }

  it should "return 400 for a malformed WDL " in {
    Post(s"/workflows/$version", FormData(Seq("wdlSource" -> CromwellApiServiceSpec.MalformedWdl, "workflowInputs" -> HelloWorld.rawInputs.toJson.toString()))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }
}
