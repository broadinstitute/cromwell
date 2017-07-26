package cromwell.backend.impl.jes

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendConfigurationDescriptor

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
      |filesystems {
      |  gcs {
      |    auth = "application-default"
      |  }
      |}
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
       |backend {
       |  allowed = ["JES"]
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
  val JesBackendConfigurationDescriptor = BackendConfigurationDescriptor(JesBackendConfig, JesGlobalConfig)
  val NoDefaultsConfigurationDescriptor = BackendConfigurationDescriptor(JesBackendNoDefaultConfig, JesGlobalConfig)
}
