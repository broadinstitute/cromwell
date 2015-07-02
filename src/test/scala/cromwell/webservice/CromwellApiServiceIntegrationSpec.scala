package cromwell.webservice

import cromwell.engine.db.DummyDataAccess
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.{FlatSpec, Matchers}
import spray.http.{FormData, StatusCodes}
import spray.json.DefaultJsonProtocol._
import spray.json._
import spray.testkit.ScalatestRouteTest


class CromwellApiServiceIntegrationSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {
  def actorRefFactory = system
  val workflowManager = system.actorOf(WorkflowManagerActor.props(DummyDataAccess()))

  it should "return 400 for a malformed WDL " in {
    Post("/workflows", FormData(Seq("wdlSource" -> CromwellApiServiceSpec.MalformedWdl , "workflowInputs" -> HelloWorld.rawInputs.toJson.toString() ))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
      }
  }
}
