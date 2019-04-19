package cromwell.backend.google.pipelines.common

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.api.PipelinesApiFactoryInterface
import cromwell.backend.google.pipelines.common.authentication.PipelinesApiDockerCredentials
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.BackendDockerConfiguration
import net.ceedubs.ficus.Ficus._
import spray.json._

import scala.concurrent.duration.FiniteDuration


class PipelinesApiConfiguration(val configurationDescriptor: BackendConfigurationDescriptor,
                                val genomicsFactory: PipelinesApiFactoryInterface,
                                val googleConfig: GoogleConfiguration,
                                val jesAttributes: PipelinesApiAttributes) extends DefaultJsonProtocol {

  val jesAuths = jesAttributes.auths
  val root = configurationDescriptor.backendConfig.getString("root")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val jesComputeServiceAccount = jesAttributes.computeServiceAccount

  val dockerCredentials = {
    BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map { creds =>
      PipelinesApiDockerCredentials.apply(creds, googleConfig)
    }
  }

  val dockerEncryptionKeyName: Option[String] = dockerCredentials flatMap { _.keyName }
  val dockerEncryptionAuthName: Option[String] = dockerCredentials flatMap { _.authName }
  val dockerToken: Option[String] = dockerCredentials map { _.token }

  val needAuthFileUpload = jesAuths.gcs.requiresAuthFile || dockerCredentials.isDefined || jesAttributes.restrictMetadataAccess
  val jobShell = configurationDescriptor.backendConfig.as[Option[String]]("job-shell").getOrElse(
    configurationDescriptor.globalConfig.getOrElse("system.job-shell", "/bin/bash"))

  lazy val slowJobWarningAfter = configurationDescriptor.backendConfig.as[Option[FiniteDuration]](path="slow-job-warning-time")
}
