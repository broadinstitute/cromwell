package cromwell.backend.google.batch.models

import akka.actor.ActorSystem
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.WorkflowOptions
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.path.PathBuilder

import scala.concurrent.Await
import scala.concurrent.duration._

object GcpBatchTestConfig {

  private val BatchBackendConfigString =
    """
      |project = "my-cromwell-workflows"
      |root = "gs://my-cromwell-workflows-bucket"
      |
      |genomics {
      |  auth = "application-default"
      |  location = "us-central1"
      |}
      |
      |filesystems.gcs.auth = "application-default"
      |
      |request-workers = 1
      |
      |default-runtime-attributes {
      |    cpu: 1
      |    failOnStderr: false
      |    continueOnReturnCode: 0
      |    docker: "ubuntu:latest"
      |    memory: "2048 MB"
      |    bootDiskSizeGb: 10
      |    disks: "local-disk 10 SSD"
      |    noAddress: false
      |    preemptible: 0
      |    zones:["us-central1-b", "us-central1-a"]
      |}
      |
      |""".stripMargin

  private val NoDefaultsConfigString =
    """
      |project = "my-cromwell-workflows"
      |root = "gs://my-cromwell-workflows-bucket"
      |
      |genomics {
      |  auth = "application-default"
      |}
      |
      |filesystems {
      |  gcs {
      |    auth = "application-default"
      |  }
      |}
      |""".stripMargin

  private val BatchGlobalConfigString =
    s"""
       |google {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = mock
       |      scheme = mock
       |    }
       |    {
       |      # legacy `application-default` auth that actually just mocks
       |      name = "application-default"
       |      scheme = "mock"
       |    }
       |  ]
       |}
       |
       |filesystems {
       |  gcs {
       |    class = "cromwell.filesystems.gcs.GcsPathBuilderFactory"
       |  }
       |}
       |
       |backend {
       |  default = "batch"
       |  providers {
       |    batch {
       |      actor-factory = "cromwell.backend.google.pipelines.batch.GcpBatchBackendLifecycleActorFactory"
       |      config {
       |      $BatchBackendConfigString
       |      }
       |    }
       |  }
       |}
       |
       |""".stripMargin

  val BatchBackendConfig: Config = ConfigFactory.parseString(BatchBackendConfigString)
  val BatchGlobalConfig: Config = ConfigFactory.parseString(BatchGlobalConfigString)
  val BatchBackendNoDefaultConfig: Config = ConfigFactory.parseString(NoDefaultsConfigString)
  val BatchBackendConfigurationDescriptor: BackendConfigurationDescriptor =
    new BackendConfigurationDescriptor(BatchBackendConfig, BatchGlobalConfig) {
      override private[backend] lazy val cromwellFileSystems = new CromwellFileSystems(BatchGlobalConfig)
    }
  val NoDefaultsConfigurationDescriptor: BackendConfigurationDescriptor =
    BackendConfigurationDescriptor(BatchBackendNoDefaultConfig, BatchGlobalConfig)
  def pathBuilders()(implicit as: ActorSystem): List[PathBuilder] =
    Await.result(BatchBackendConfigurationDescriptor.pathBuilders(WorkflowOptions.empty), 5.seconds)
  val googleConfiguration: GoogleConfiguration = GoogleConfiguration(BatchGlobalConfig)
  val batchAttributes: GcpBatchConfigurationAttributes =
    GcpBatchConfigurationAttributes(googleConfiguration, BatchBackendConfig, "batch")
  val gcpBatchConfiguration =
    new GcpBatchConfiguration(BatchBackendConfigurationDescriptor, googleConfiguration, batchAttributes)
}
