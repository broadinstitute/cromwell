package cromwell.backend.google.pipelines.v2beta.api.request

import java.net.URL
import java.time.OffsetDateTime
import akka.actor.ActorRef
import com.google.api.client.http.GenericUrl
import com.google.api.client.testing.http.MockHttpTransport
import com.google.api.services.lifesciences.v2beta.model.Operation
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIStatusPollRequest
import cromwell.backend.google.pipelines.common.api.RunStatus._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.{ExecutionEvent, WorkflowId}
import io.grpc.Status
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class GetRequestHandlerSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "GetRequestHandler"

  private val requestHandler: GetRequestHandler = new RequestHandler(
    "GetRequestHandlerSpec",
    new URL("file:///getrequesthandlerspec"),
    BatchRequestTimeoutConfiguration(None, None),
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
    ("check that we classify error code 10 as a preemption on a preemptible VM",
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
        |}""".stripMargin,
      Preempted(
        Status.ABORTED,
        None,
        Nil,
        List(
          ExecutionEvent("waiting for quota", OffsetDateTime.parse("2019-08-18T12:04:38.082650Z"),None),
          ExecutionEvent("Complete in GCE / Cromwell Poll Interval", OffsetDateTime.parse("2019-08-18T15:58:26.659602622Z"),None),
        ),
        Some("custom-2-7168"),
        None,
        None)
    ),
    ("check that we classify error code 10 as a failure on a non-preemptible VM",
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
        |          "preemptible": false,
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
        |}""".stripMargin,
      Failed(
        Status.ABORTED,
        None,
        Nil,
        List(
          ExecutionEvent("waiting for quota", OffsetDateTime.parse("2019-08-18T12:04:38.082650Z"),None),
          ExecutionEvent("Complete in GCE / Cromwell Poll Interval", OffsetDateTime.parse("2019-08-18T15:58:26.659602622Z"),None),
        ),
        Some("custom-2-7168"),
        None,
        None
      )
    ),
    // As of 2022-01 the zone `us-west3` in `broad-dsde-cromwell-dev` has its CPU quota purposely de-rated to 1 for testing
    ("check that a job is AwaitingCloudQuota if its most recent event is quota exhaustion",
      """{
        |  "metadata": {
        |    "@type": "type.googleapis.com/google.cloud.lifesciences.v2beta.Metadata",
        |    "createTime": "2022-01-19T21:53:55.138960Z",
        |    "events": [
        |      {
        |        "delayed": {
        |          "cause": "generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high",
        |          "metrics": [
        |            "CPUS"
        |          ]
        |        },
        |        "description": "A resource limit has delayed the operation: generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high",
        |        "timestamp": "2022-01-19T21:54:07.717679160Z"
        |      }
        |    ],
        |    "labels": {
        |      "cromwell-workflow-id": "cromwell-ac888b4e-2e6b-4dcc-a537-3c6db7764037",
        |      "wdl-task-name": "sleep"
        |    },
        |    "pipeline": {
        |      "actions": [
        |        {
        |          "commands": [
        |            "-c",
        |            "printf '%s %s\\n' \"$(date -u '+%Y/%m/%d %H:%M:%S')\" Starting\\ container\\ setup."
        |          ],
        |          "entrypoint": "/bin/sh",
        |          "imageUri": "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine",
        |          "labels": {
        |            "logging": "ContainerSetup"
        |          },
        |          "timeout": "300s"
        |        }
        |      ],
        |      "environment": {
        |        "MEM_SIZE": "2.0",
        |        "MEM_UNIT": "GB"
        |      },
        |      "resources": {
        |        "virtualMachine": {
        |          "bootDiskSizeGb": 12,
        |          "bootImage": "projects/cos-cloud/global/images/family/cos-stable",
        |          "disks": [
        |            {
        |              "name": "local-disk",
        |              "sizeGb": 10,
        |              "type": "pd-ssd"
        |            }
        |          ],
        |          "labels": {
        |            "cromwell-workflow-id": "cromwell-ac888b4e-2e6b-4dcc-a537-3c6db7764037",
        |            "goog-pipelines-worker": "true",
        |            "wdl-task-name": "sleep"
        |          },
        |          "machineType": "custom-1-2048",
        |          "network": {},
        |          "nvidiaDriverVersion": "450.51.06",
        |          "serviceAccount": {
        |            "email": "centaur@broad-dsde-cromwell-dev.iam.gserviceaccount.com",
        |            "scopes": [
        |              "https://www.googleapis.com/auth/compute",
        |              "https://www.googleapis.com/auth/devstorage.full_control",
        |              "https://www.googleapis.com/auth/cloudkms",
        |              "https://www.googleapis.com/auth/userinfo.email",
        |              "https://www.googleapis.com/auth/userinfo.profile",
        |              "https://www.googleapis.com/auth/monitoring.write",
        |              "https://www.googleapis.com/auth/bigquery",
        |              "https://www.googleapis.com/auth/cloud-platform"
        |            ]
        |          },
        |          "volumes": [
        |            {
        |              "persistentDisk": {
        |                "sizeGb": 10,
        |                "type": "pd-ssd"
        |              },
        |              "volume": "local-disk"
        |            }
        |          ]
        |        },
        |        "zones": [
        |          "us-west3-a",
        |          "us-west3-b",
        |          "us-west3-c"
        |        ]
        |      },
        |      "timeout": "604800s"
        |    }
        |  },
        |  "name": "projects/1005074806481/locations/us-central1/operations/3874882033889365536"
        |}""".stripMargin,
      AwaitingCloudQuota
    ),
    ("check that a job is Running and no longer AwaitingCloudQuota once a worker assigns",
      """{
        |  "metadata": {
        |    "@type": "type.googleapis.com/google.cloud.lifesciences.v2beta.Metadata",
        |    "createTime": "2022-01-19T21:53:55.138960Z",
        |    "events": [
        |      {
        |        "description": "Started pulling \"gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine\"",
        |        "pullStarted": {
        |          "imageUri": "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine"
        |        },
        |        "timestamp": "2022-01-19T22:09:55.410251187Z"
        |      },
        |      {
        |        "description": "Worker \"google-pipelines-worker-e6c8bf8035860b2cd69488497bd602d8\" assigned in \"us-west3-c\" on a \"custom-1-2048\" machine",
        |        "timestamp": "2022-01-19T22:09:20.363771714Z",
        |        "workerAssigned": {
        |          "instance": "google-pipelines-worker-e6c8bf8035860b2cd69488497bd602d8",
        |          "machineType": "custom-1-2048",
        |          "zone": "us-west3-c"
        |        }
        |      },
        |      {
        |        "delayed": {
        |          "cause": "generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high",
        |          "metrics": [
        |            "CPUS"
        |          ]
        |        },
        |        "description": "A resource limit has delayed the operation: generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high",
        |        "timestamp": "2022-01-19T21:54:07.717679160Z"
        |      }
        |    ],
        |    "labels": {
        |      "cromwell-workflow-id": "cromwell-ac888b4e-2e6b-4dcc-a537-3c6db7764037",
        |      "wdl-task-name": "sleep"
        |    },
        |    "pipeline": {
        |      "actions": [
        |        {
        |          "commands": [
        |            "-c",
        |            "printf '%s %s\\n' \"$(date -u '+%Y/%m/%d %H:%M:%S')\" Starting\\ container\\ setup."
        |          ],
        |          "entrypoint": "/bin/sh",
        |          "imageUri": "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine",
        |          "labels": {
        |            "logging": "ContainerSetup"
        |          },
        |          "timeout": "300s"
        |        }
        |      ],
        |      "environment": {
        |        "MEM_SIZE": "2.0",
        |        "MEM_UNIT": "GB"
        |      },
        |      "resources": {
        |        "virtualMachine": {
        |          "bootDiskSizeGb": 12,
        |          "bootImage": "projects/cos-cloud/global/images/family/cos-stable",
        |          "disks": [
        |            {
        |              "name": "local-disk",
        |              "sizeGb": 10,
        |              "type": "pd-ssd"
        |            }
        |          ],
        |          "labels": {
        |            "cromwell-workflow-id": "cromwell-ac888b4e-2e6b-4dcc-a537-3c6db7764037",
        |            "goog-pipelines-worker": "true",
        |            "wdl-task-name": "sleep"
        |          },
        |          "machineType": "custom-1-2048",
        |          "network": {},
        |          "nvidiaDriverVersion": "450.51.06",
        |          "serviceAccount": {
        |            "email": "centaur@broad-dsde-cromwell-dev.iam.gserviceaccount.com",
        |            "scopes": [
        |              "https://www.googleapis.com/auth/compute",
        |              "https://www.googleapis.com/auth/devstorage.full_control",
        |              "https://www.googleapis.com/auth/cloudkms",
        |              "https://www.googleapis.com/auth/userinfo.email",
        |              "https://www.googleapis.com/auth/userinfo.profile",
        |              "https://www.googleapis.com/auth/monitoring.write",
        |              "https://www.googleapis.com/auth/bigquery",
        |              "https://www.googleapis.com/auth/cloud-platform"
        |            ]
        |          },
        |          "volumes": [
        |            {
        |              "persistentDisk": {
        |                "sizeGb": 10,
        |                "type": "pd-ssd"
        |              },
        |              "volume": "local-disk"
        |            }
        |          ]
        |        },
        |        "zones": [
        |          "us-west3-a",
        |          "us-west3-b",
        |          "us-west3-c"
        |        ]
        |      },
        |      "timeout": "604800s"
        |    },
        |    "startTime": "2022-01-19T22:09:20.363771714Z"
        |  },
        |  "name": "projects/1005074806481/locations/us-central1/operations/3874882033889365536"
        |}
        |
        |
        |""".stripMargin,
      Running
    ),
    ("check that a job is no longer AwaitingCloudQuota once it finishes",
      """{
        |  "done": true,
        |  "metadata": {
        |    "@type": "type.googleapis.com/google.cloud.lifesciences.v2beta.Metadata",
        |    "createTime": "2022-01-19T19:17:13.175579Z",
        |    "endTime": "2022-01-19T19:37:22.764120036Z",
        |    "events": [
        |      {
        |        "description": "Worker released",
        |        "timestamp": "2022-01-19T19:37:22.764120036Z",
        |        "workerReleased": {
        |          "instance": "google-pipelines-worker-8eff543e6858c204c8f67520aee75432",
        |          "zone": "us-west3-c"
        |        }
        |      },
        |      {
        |        "containerStopped": {
        |          "actionId": 19
        |        },
        |        "description": "Stopped running shortened for test",
        |        "timestamp": "2022-01-19T19:37:19.822873814Z"
        |      },
        |      {
        |        "description": "Started pulling \"gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine\"",
        |        "pullStarted": {
        |          "imageUri": "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine"
        |        },
        |        "timestamp": "2022-01-19T19:32:55.709674372Z"
        |      },
        |      {
        |        "description": "Worker \"google-pipelines-worker-8eff543e6858c204c8f67520aee75432\" assigned in \"us-west3-c\" on a \"custom-1-2048\" machine",
        |        "timestamp": "2022-01-19T19:32:19.204055448Z",
        |        "workerAssigned": {
        |          "instance": "google-pipelines-worker-8eff543e6858c204c8f67520aee75432",
        |          "machineType": "custom-1-2048",
        |          "zone": "us-west3-c"
        |        }
        |      },
        |      {
        |        "delayed": {
        |          "cause": "generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high",
        |          "metrics": [
        |            "CPUS"
        |          ]
        |        },
        |        "description": "A resource limit has delayed the operation: generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high",
        |        "timestamp": "2022-01-19T19:17:14.948193837Z"
        |      }
        |    ],
        |    "labels": {
        |      "cromwell-workflow-id": "cromwell-058bff35-4a55-4c0f-9113-0885f4119cd9",
        |      "wdl-task-name": "sleep"
        |    },
        |    "pipeline": {
        |      "actions": [
        |        {
        |          "commands": [
        |            "-c",
        |            "printf '%s %s\\n' \"$(date -u '+%Y/%m/%d %H:%M:%S')\" Starting\\ container\\ setup."
        |          ],
        |          "entrypoint": "/bin/sh",
        |          "imageUri": "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine",
        |          "labels": {
        |            "logging": "ContainerSetup"
        |          },
        |          "timeout": "300s"
        |        },
        |        {
        |          "alwaysRun": true,
        |          "commands": [
        |            "-c",
        |            "python3 -c 'import base64; shortened for test"
        |          ],
        |          "entrypoint": "/bin/sh",
        |          "imageUri": "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine",
        |          "labels": {
        |            "tag": "Delocalization"
        |          }
        |        }
        |      ],
        |      "environment": {
        |        "MEM_SIZE": "2.0",
        |        "MEM_UNIT": "GB"
        |      },
        |      "resources": {
        |        "virtualMachine": {
        |          "bootDiskSizeGb": 12,
        |          "bootImage": "projects/cos-cloud/global/images/family/cos-stable",
        |          "disks": [
        |            {
        |              "name": "local-disk",
        |              "sizeGb": 10,
        |              "type": "pd-ssd"
        |            }
        |          ],
        |          "labels": {
        |            "cromwell-workflow-id": "cromwell-058bff35-4a55-4c0f-9113-0885f4119cd9",
        |            "goog-pipelines-worker": "true",
        |            "wdl-task-name": "sleep"
        |          },
        |          "machineType": "custom-1-2048",
        |          "network": {},
        |          "nvidiaDriverVersion": "450.51.06",
        |          "serviceAccount": {
        |            "email": "centaur@broad-dsde-cromwell-dev.iam.gserviceaccount.com",
        |            "scopes": [
        |              "https://www.googleapis.com/auth/compute",
        |              "https://www.googleapis.com/auth/devstorage.full_control",
        |              "https://www.googleapis.com/auth/cloudkms",
        |              "https://www.googleapis.com/auth/userinfo.email",
        |              "https://www.googleapis.com/auth/userinfo.profile",
        |              "https://www.googleapis.com/auth/monitoring.write",
        |              "https://www.googleapis.com/auth/bigquery",
        |              "https://www.googleapis.com/auth/cloud-platform"
        |            ]
        |          },
        |          "volumes": [
        |            {
        |              "persistentDisk": {
        |                "sizeGb": 10,
        |                "type": "pd-ssd"
        |              },
        |              "volume": "local-disk"
        |            }
        |          ]
        |        },
        |        "zones": [
        |          "us-west3-a",
        |          "us-west3-b",
        |          "us-west3-c"
        |        ]
        |      },
        |      "timeout": "604800s"
        |    },
        |    "startTime": "2022-01-19T19:32:19.204055448Z"
        |  },
        |  "name": "projects/1005074806481/locations/us-central1/operations/5001350794958839237",
        |  "response": {
        |    "@type": "type.googleapis.com/cloud.lifesciences.pipelines.RunPipelineResponse"
        |  }
        |}""".stripMargin,
      Success(List(
        new ExecutionEvent("waiting for quota", OffsetDateTime.parse("2022-01-19T19:17:13.175579Z"), None),
        new ExecutionEvent("Worker released", OffsetDateTime.parse("2022-01-19T19:37:22.764120036Z"), None),
        new ExecutionEvent("Stopped running shortened for test", OffsetDateTime.parse("2022-01-19T19:37:19.822873814Z"), None),
        new ExecutionEvent("Started pulling \"gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine\"", OffsetDateTime.parse("2022-01-19T19:32:55.709674372Z"), Option("Pulling \"gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine\"")),
        new ExecutionEvent("Worker \"google-pipelines-worker-8eff543e6858c204c8f67520aee75432\" assigned in \"us-west3-c\" on a \"custom-1-2048\" machine", OffsetDateTime.parse("2022-01-19T19:32:19.204055448Z"), None),
        new ExecutionEvent("A resource limit has delayed the operation: generic::resource_exhausted: allocating: selecting resources: selecting region and zone: no available zones: us-west3: 1 CPUS (0/1 available) usage too high", OffsetDateTime.parse("2022-01-19T19:17:14.948193837Z"), None),
        new ExecutionEvent("Complete in GCE / Cromwell Poll Interval", OffsetDateTime.parse("2022-01-19T19:37:22.764120036Z"), None)
      ), Option("custom-1-2048"), Option("us-west3-c"), Option("google-pipelines-worker-8eff543e6858c204c8f67520aee75432"))
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
