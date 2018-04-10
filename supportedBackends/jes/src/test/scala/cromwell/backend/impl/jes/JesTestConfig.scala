package cromwell.backend.impl.jes

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.core.WorkflowOptions
import cromwell.core.filesystem.CromwellFileSystems

import scala.concurrent.Await
import scala.concurrent.duration._
object JesTestConfig {

  private val JesBackendConfigString =
    """
      |project = "my-cromwell-workflows"
      |root = "gs://my-cromwell-workflows-bucket"
      |
      |genomics {
      |  auth = "application-default"
      |  endpoint-url = "https://genomics.googleapis.com/"
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
      |    memory: "2 GB"
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
      |  endpoint-url = "https://genomics.googleapis.com/"
      |}
      |
      |filesystems {
      |  gcs {
      |    auth = "application-default"
      |  }
      |}
      |""".stripMargin

  private val JesGlobalConfigString =
    s"""
       |google {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = "application-default"
       |      scheme = "application_default"
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
       |  default = "JES"
       |  providers {
       |    JES {
       |      actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleFactory"
       |      config {
       |      $JesBackendConfigString
       |      }
       |    }
       |  }
       |}
       |
       |""".stripMargin

  val JesBackendConfig = ConfigFactory.parseString(JesBackendConfigString)
  val JesGlobalConfig = ConfigFactory.parseString(JesGlobalConfigString)
  val JesBackendNoDefaultConfig = ConfigFactory.parseString(NoDefaultsConfigString)
  val mockFilesystems = new CromwellFileSystems(JesGlobalConfig)
  val JesBackendConfigurationDescriptor = new BackendConfigurationDescriptor(JesBackendConfig, JesGlobalConfig) {
    override lazy val configuredPathBuilderFactories = mockFilesystems.factoriesFromConfig(JesBackendConfig).unsafe("Failed to instantiate backend filesystem")
  }
  def pathBuilders()(implicit as: ActorSystem) = Await.result(JesBackendConfigurationDescriptor.pathBuilders(WorkflowOptions.empty), 5.seconds)
  val NoDefaultsConfigurationDescriptor = BackendConfigurationDescriptor(JesBackendNoDefaultConfig, JesGlobalConfig)
}
