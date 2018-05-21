package cromwell.backend.google.pipelines.common

import akka.actor.ActorSystem
import com.google.api.client.http.HttpRequestInitializer
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.WorkflowOptions
import cromwell.core.filesystem.CromwellFileSystems

import scala.concurrent.Await
import scala.concurrent.duration._

object PipelinesApiTestConfig {

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
       |      actor-factory = "cromwell.backend.google.pipelines.common.PipelinesApiBackendLifecycleActorFactory"
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
  val JesBackendConfigurationDescriptor = new BackendConfigurationDescriptor(JesBackendConfig, JesGlobalConfig) {
    override private[backend] lazy val cromwellFileSystems = new CromwellFileSystems(JesGlobalConfig)
  }
  val NoDefaultsConfigurationDescriptor = BackendConfigurationDescriptor(JesBackendNoDefaultConfig, JesGlobalConfig)
  val genomicsFactory = new PipelinesApiFactoryInterface {
    override def build(httpRequestInitializer: HttpRequestInitializer) = new PipelinesApiRequestFactory {
      override def cancelRequest(job: StandardAsyncJob) = ???
      override def getRequest(job: StandardAsyncJob) = ???
      override def runRequest(createPipelineParameters: PipelinesApiRequestFactory.CreatePipelineParameters, jobPath: Option[JobPaths] = None) = ???
    }
  }
  def pathBuilders()(implicit as: ActorSystem) = Await.result(JesBackendConfigurationDescriptor.pathBuilders(WorkflowOptions.empty), 5.seconds)
  val googleConfiguration = GoogleConfiguration(JesGlobalConfig)
  val jesAttributes = PipelinesApiAttributes(googleConfiguration, JesBackendConfig)
  val jesConfiguration = new PipelinesApiConfiguration(JesBackendConfigurationDescriptor, genomicsFactory, googleConfiguration, jesAttributes)
}
