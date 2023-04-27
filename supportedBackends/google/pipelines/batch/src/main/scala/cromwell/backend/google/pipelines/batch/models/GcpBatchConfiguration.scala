package cromwell.backend.google.pipelines.batch.models

import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.batch.authentication.{GcpBatchAuths, GcpBatchDockerCredentials}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.BackendDockerConfiguration
import net.ceedubs.ficus.Ficus._
import spray.json._

import scala.concurrent.duration.FiniteDuration


class GcpBatchConfiguration(val configurationDescriptor: BackendConfigurationDescriptor,
                            val googleConfig: GoogleConfiguration,
                            val batchAttributes: GcpBatchConfigurationAttributes
                           ) extends DefaultJsonProtocol {

  val batchAuths: GcpBatchAuths = batchAttributes.auths
  val root: String = configurationDescriptor.backendConfig.getString("root")
  val batchTimeout: FiniteDuration = batchAttributes.batchTimeout
  val runtimeConfig: Option[Config] = configurationDescriptor.backendRuntimeAttributesConfig


  val dockerCredentials: Option[GcpBatchDockerCredentials] = {
    BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map { creds =>
      GcpBatchDockerCredentials.apply(creds, googleConfig)
   }
  }

  val dockerEncryptionKeyName: Option[String] = dockerCredentials flatMap { _.keyName }
  val dockerEncryptionAuthName: Option[String] = dockerCredentials flatMap { _.authName }
  val dockerToken: Option[String] = dockerCredentials map { _.token }

  val jobShell: String = configurationDescriptor.backendConfig.as[Option[String]]("job-shell").getOrElse(
    configurationDescriptor.globalConfig.getOrElse("system.job-shell", "/bin/bash"))

}
