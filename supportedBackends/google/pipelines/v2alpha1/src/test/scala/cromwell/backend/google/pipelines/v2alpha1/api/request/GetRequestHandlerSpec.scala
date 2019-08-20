package cromwell.backend.google.pipelines.v2alpha1.api.request

import java.net.URL
import java.time.OffsetDateTime

import akka.actor.ActorRef
import com.google.api.client.http.GenericUrl
import com.google.api.client.testing.http.MockHttpTransport
import com.google.api.services.genomics.v2alpha1.model.Operation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIStatusPollRequest
import cromwell.backend.google.pipelines.common.api.RunStatus._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.{ExecutionEvent, WorkflowId}
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
    ("Check that we classify error code 10 as a preemption",
      """{
        |  "done": true,
        |  "error": {
        |    "code": 10,
        |    "message": "The assigned worker has failed to complete the operation"
        |  },
        |  "metadata": {
        |    "@type": "type.googleapis.com/google.genomics.v2alpha1.Metadata",
        |    "createTime": "2019-08-18T12:04:38.082650Z",
        |    "endTime": "2019-08-18T15:58:26.659602622Z",
        |    "events": [],
        |    "labels": {
        |      "cromwell-sub-workflow-name": "bamtocram",
        |      "cromwell-workflow-id": "asdfasdf",
        |      "wdl-call-alias": "validatecram",
        |      "wdl-task-name": "validatesamfile"
        |    },
        |    "pipeline": {
        |      "actions": [],
        |      "environment": {},
        |      "resources": {
        |        "projectId": "",
        |        "regions": [],
        |        "virtualMachine": {
        |          "accelerators": [],
        |          "bootDiskSizeGb": 11,
        |          "bootImage": "asdfasdf",
        |          "cpuPlatform": "",
        |          "disks": [
        |            {
        |              "name": "local-disk",
        |              "sizeGb": 41,
        |              "sourceImage": "",
        |              "type": "pd-standard"
        |            }
        |          ],
        |          "enableStackdriverMonitoring": false,
        |          "labels": {
        |            "cromwell-sub-workflow-name": "bamtocram",
        |            "cromwell-workflow-id": "asdfasdf",
        |            "goog-pipelines-worker": "true",
        |            "wdl-call-alias": "validatecram",
        |            "wdl-task-name": "validatesamfile"
        |          },
        |          "machineType": "custom-2-7168",
        |          "network": {
        |            "name": "",
        |            "subnetwork": "",
        |            "usePrivateAddress": false
        |          },
        |          "nvidiaDriverVersion": "",
        |          "preemptible": true,
        |          "serviceAccount": {
        |            "email": "default",
        |            "scopes": [
        |              "https://www.googleapis.com/auth/genomics",
        |              "https://www.googleapis.com/auth/compute",
        |              "https://www.googleapis.com/auth/devstorage.full_control",
        |              "https://www.googleapis.com/auth/cloudkms",
        |              "https://www.googleapis.com/auth/userinfo.email",
        |              "https://www.googleapis.com/auth/userinfo.profile",
        |              "https://www.googleapis.com/auth/monitoring.write",
        |              "https://www.googleapis.com/auth/cloud-platform"
        |            ]
        |          }
        |        },
        |        "zones": [
        |          "us-central1-a",
        |          "us-central1-b",
        |          "us-east1-d",
        |          "us-central1-c",
        |          "us-central1-f",
        |          "us-east1-c"
        |        ]
        |      },
        |      "timeout": "604800s"
        |    },
        |    "startTime": "2019-08-18T12:04:39.192909594Z"
        |  },
        |  "name": "asdfasdf"
        |}""".stripMargin, Preempted(Status.ABORTED, None, Nil, List(ExecutionEvent("waiting for quota", OffsetDateTime.parse("2019-08-18T12:04:38.082650Z"),None)), Some("custom-2-7168"), None, None)
    )
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
