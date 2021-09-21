package cromwell.backend.google.pipelines.common

import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.api.PipelinesApiFactoryInterface
import cromwell.backend.google.pipelines.common.authentication.{PipelinesApiAuths, PipelinesApiDockerCredentials}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.BackendDockerConfiguration
import net.ceedubs.ficus.Ficus._
import spray.json._

import scala.concurrent.duration.FiniteDuration

class PipelinesApiConfiguration(val configurationDescriptor: BackendConfigurationDescriptor,
                                val genomicsFactory: PipelinesApiFactoryInterface,
                                val googleConfig: GoogleConfiguration,
                                val papiAttributes: PipelinesApiConfigurationAttributes) extends DefaultJsonProtocol {

  val jesAuths: PipelinesApiAuths = papiAttributes.auths
  val root: String = configurationDescriptor.backendConfig.getString("root")
  val pipelineTimeout: FiniteDuration = papiAttributes.pipelineTimeout
  val runtimeConfig: Option[Config] = configurationDescriptor.backendRuntimeAttributesConfig

  val dockerCredentials: Option[PipelinesApiDockerCredentials] = {
    BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map { creds =>
      PipelinesApiDockerCredentials.apply(creds, googleConfig)
    }
  }

  val dockerEncryptionKeyName: Option[String] = dockerCredentials flatMap { _.keyName }
  val dockerEncryptionAuthName: Option[String] = dockerCredentials flatMap { _.authName }
  val dockerToken: Option[String] = dockerCredentials map { _.token }

  val jobShell: String = configurationDescriptor.backendConfig.as[Option[String]]("job-shell").getOrElse(
    configurationDescriptor.globalConfig.getOrElse("system.job-shell", "/bin/bash"))
}
