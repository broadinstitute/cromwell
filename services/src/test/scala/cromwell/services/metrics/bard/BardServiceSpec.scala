package cromwell.services.metrics.bard

import cromwell.services.metrics.bard.model.TaskSummaryEvent
import org.mockserver.integration.ClientAndServer
import org.mockserver.model.{HttpRequest, HttpResponse, JsonBody, MediaType}
import org.mockserver.verify.VerificationTimes
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class BardServiceSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "sendEvent in BardService"

  it should "reach out to Bard to send an event" in {
    val mockServer = ClientAndServer.startClientAndServer()
    try {
      val bardApiPath = "/api/eventLog/v1/cromwell/task:summary"
      val bardUrl = "http://localhost:" + mockServer.getPort
      //  Mock the server response
      mockServer
        .when(HttpRequest.request(bardApiPath).withMethod("POST"))
        .respond(HttpResponse.response.withStatusCode(200).withContentType(MediaType.APPLICATION_JSON))

      val bardService = new BardService(bardUrl, 10)
      bardService.sendEvent("bearerToken", TaskSummaryEvent("done", 0))

      mockServer.verify(
        HttpRequest
          .request(bardApiPath)
          .withMethod("POST")
          .withBody(new JsonBody("""{
                                   |   "properties": {
                                   |     "terminationStatus": "done",
                                   |     "returnCode": 0,
                                   |     "event": "task:summary",
                                   |     "pushToMixpanel": false
                                   |   }
                                   | }
                                   |""".stripMargin)),
        VerificationTimes.exactly(1)
      )
    } finally if (mockServer != null) mockServer.close()
  }

  it should "not throw if Bard returns a failed response" in {
    val mockServer = ClientAndServer.startClientAndServer()
    try {
      val bardApiPath = "/api/eventLog/v1/cromwell/task:summary"
      val bardUrl = "http://localhost:" + mockServer.getPort
      //  Mock the server response
      mockServer
        .when(HttpRequest.request(bardApiPath).withMethod("POST"))
        .respond(HttpResponse.response.withStatusCode(500).withContentType(MediaType.APPLICATION_JSON))

      val bardService = new BardService(bardUrl, 10)
      bardService.sendEvent("bearerToken", TaskSummaryEvent("done", 0))
    } finally if (mockServer != null) mockServer.close()
  }

}
