package cromwell.services.metrics.bard

import org.mockserver.integration.ClientAndServer
import org.mockserver.model.{HttpRequest, HttpResponse, JsonBody, MediaType}
import org.mockserver.verify.VerificationTimes
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class BardServiceSpec extends AnyFlatSpec with Matchers with BardTestUtils {

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

      val bardService = new BardService(bardUrl, 10, serviceRegistryProbe.ref)

      bardService.sendEvent(taskSummaryEvent)

      mockServer.verify(
        HttpRequest
          .request(bardApiPath)
          .withMethod("POST")
          .withBody(new JsonBody(s"""{
                                    |   "properties": {
                                    |     "workflowId": "$workflowId",
                                    |     "rootWorkflowId": "$rootWorkflowId",
                                    |     "jobTag": "$jobTag",
                                    |     "jobFullyQualifiedName": "$jobFqn",
                                    |     "jobIndex": $jobIndex,
                                    |     "jobAttempt": $jobAttempt,
                                    |     "terminalState": "$terminalState",
                                    |     "platform": "$platform",
                                    |     "dockerImage": "$dockerImage",
                                    |     "cpuCount": $cpu,
                                    |     "memoryBytes": $memory,
                                    |     "startTime": "$start",
                                    |     "cpuStartTime": "$cpuStart",
                                    |     "endTime": "$end",
                                    |     "jobSeconds": $jobSeconds,
                                    |     "cpuSeconds": $cpuSeconds,
                                    |     "distinct_id": "$workflowId"
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
        .respond(
          HttpResponse.response
            .withStatusCode(500)
            .withContentType(MediaType.APPLICATION_JSON)
            .withBody("""{ "error": "Expected Test Error" }""")
        )

      val bardService = new BardService(bardUrl, 10, serviceRegistryProbe.ref)

      bardService.sendEvent(taskSummaryEvent)

    } finally if (mockServer != null) mockServer.close()
  }

}
