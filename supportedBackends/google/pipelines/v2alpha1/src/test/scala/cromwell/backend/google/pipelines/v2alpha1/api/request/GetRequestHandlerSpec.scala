package cromwell.backend.google.pipelines.v2alpha1.api.request

import java.net.URL

import akka.actor.ActorRef
import com.google.api.client.http.GenericUrl
import com.google.api.client.testing.http.MockHttpTransport
import com.google.api.services.genomics.v2alpha1.model.Operation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIStatusPollRequest
import cromwell.backend.google.pipelines.common.api.RunStatus._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.WorkflowId
import io.grpc.Status
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

class GetRequestHandlerSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "GetRequestHandler"

  private val requestHandler: GetRequestHandler = new RequestHandler(
    "GetRequestHandlerSpec",
    new URL("file:///getrequesthandlerspec"),
    null // This can be null because we never need to initialize http requests for this test
  )

  private val workflowId = WorkflowId.randomId()
  private val httpRequest =
    new MockHttpTransport.Builder().build().createRequestFactory().buildGetRequest(new GenericUrl())
  private val actorRef = ActorRef.noSender
  private val jobId = StandardAsyncJob("test_job")
  private val pollingRequest = PAPIStatusPollRequest(workflowId, actorRef, httpRequest, jobId)

  private val interpretedStatus = Table(
    ("description", "json", "status"),
    ("parse null operation json", null, UnsuccessfulRunStatus(
      Status.UNKNOWN,
      Option("Operation returned as empty"),
      Nil,
      None,
      None,
      None,
      wasPreemptible = false
    )),
    ("parse empty operation json", "{}", Initializing),
    ("parse error operation json without resources",
      """|{
         |  "done": true,
         |  "error": {}
         |}
         |""".stripMargin,
      Failed(Status.UNAVAILABLE, None, Nil, Nil, None, None, None)
    ),
    ("parse error operation json without virtualMachine",
      """|{
         |  "done": true,
         |  "resources": {
         |  },
         |  "error": {}
         |}
         |""".stripMargin,
      Failed(Status.UNAVAILABLE, None, Nil, Nil, None, None, None)
    ),
    ("parse error operation json without preemptible",
      """|{
         |  "done": true,
         |  "resources": {
         |    "virtualMachine": {
         |    }
         |  },
         |  "error": {}
         |}
         |""".stripMargin,
      Failed(Status.UNAVAILABLE, None, Nil, Nil, None, None, None)
    ),
    ("parse error operation json with preemptible true",
      """|{
         |  "done": true,
         |  "resources": {
         |    "virtualMachine": {
         |      "preemptible": true
         |    }
         |  },
         |  "error": {}
         |}
         |""".stripMargin,
      Failed(Status.UNAVAILABLE, None, Nil, Nil, None, None, None)
    ),
    ("parse error operation json with preemptible false",
      """|{
         |  "done": true,
         |  "resources": {
         |    "virtualMachine": {
         |      "preemptible": false
         |    }
         |  },
         |  "error": {}
         |}
         |""".stripMargin,
      Failed(Status.UNAVAILABLE, None, Nil, Nil, None, None, None)
    ),
  )

  forAll(interpretedStatus) { (description, json, expectedStatus) =>
    it should description in {
      // Operation responses could come back as null. Handle it and don't crash.
      // https://github.com/googleapis/google-http-java-client/blob/v1.28.0/google-http-client/src/main/java/com/google/api/client/http/HttpResponse.java#L456-L458
      val operation =
      Option(json).map(GoogleAuthMode.jsonFactory.createJsonParser).map(_.parse(classOf[Operation])).orNull
      val runStatus = requestHandler.interpretOperationStatus(operation, pollingRequest)
      runStatus should be(expectedStatus)
    }
  }

}
