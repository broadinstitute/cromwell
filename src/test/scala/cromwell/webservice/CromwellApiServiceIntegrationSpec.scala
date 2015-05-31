package cromwell.webservice

import cromwell.engine.ActorWorkflowManager
import org.scalatest.{Matchers, FlatSpec}
import spray.http.{StatusCodes, FormData}
import spray.testkit.ScalatestRouteTest


class CromwellApiServiceIntegrationSpec extends FlatSpec with CromwellApiService with ScalatestRouteTest with Matchers {

  val workflowManagerActorRef = system.actorOf(ActorWorkflowManager.props)

  def actorRefFactory = system

  it should "return 400 for a malformed WDL " in {
    Post("/workflows", FormData(Seq("wdlSource" -> CromwellApiServiceSpec.MalformedWdl , "workflowInputs" -> CromwellApiServiceSpec.HelloInputsJson ))) ~>
      submitRoute ~>
      check {
        assertResult(StatusCodes.BadRequest) {
          status
        }
        assertResult("wdlSource was malformed") {
          responseAs[String]
        }
      }
  }

  // TODO: add test to check if (a) json parsing of parameter map fails (b) required fields missing

}
