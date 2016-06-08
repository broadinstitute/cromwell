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
      |""".stripMargin

  private val JesGlobalConfigString =
    s"""
       |google {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = "application-default"
       |      scheme = "application_default"
       |    },
       |    {
       |      name = "user-via-refresh"
       |      scheme = "refresh_token"
       |      client-id = "secret_id"
       |      client-secret = "secret_secret"
       |    },
       |    {
       |      name = "service-account"
       |      scheme = "service_account"
       |      service-account-id = "my-service-account"
       |      pem-file = "/path/to/file.pem"
       |    }
       |  ]
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
       |""".stripMargin

  val JesBackendConfig = ConfigFactory.parseString(JesBackendConfigString)
  val JesGlobalConfig = ConfigFactory.parseString(JesGlobalConfigString)
  val JesBackendConfigurationDescriptor = new BackendConfigurationDescriptor(JesBackendConfig, JesGlobalConfig)
}
