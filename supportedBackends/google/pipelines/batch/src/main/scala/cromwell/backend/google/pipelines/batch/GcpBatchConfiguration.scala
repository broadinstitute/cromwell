package cromwell.backend.google.pipelines.batch

//import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.cloudsupport.gcp.GoogleConfiguration

//import scala.concurrent.duration.FiniteDuration
//import cromwell.backend.google.pipelines.batch.authentication.GcpBatchAuths
import net.ceedubs.ficus.Ficus._
import spray.json._

//import scala.concurrent.duration.FiniteDuration

class GcpBatchConfiguration(val configurationDescriptor: BackendConfigurationDescriptor,
                            //val batchFactory: GcpBatchFactoryInterface,
                            val googleConfig: GoogleConfiguration,
                            //val batchAttributes: GcpBatchConfigurationAttributes
                           ) extends DefaultJsonProtocol {

  //val batchAuths: GcpBatchAuths = batchAttributes.auths
  val root: String = configurationDescriptor.backendConfig.getString("root")
  //val pipelineTimeout: FiniteDuration = batchAttributes.pipelineTimeout
  //val runtimeConfig: Option[Config] = configurationDescriptor.backendRuntimeAttributesConfig
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig

  //val dockerCredentials: Option[PipelinesApiDockerCredentials] = {
  //  BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map { creds =>
  //    PipelinesApiDockerCredentials.apply(creds, googleConfig)
  // }
  //}

  //val dockerEncryptionKeyName: Option[String] = dockerCredentials flatMap { _.keyName }
  //val dockerEncryptionAuthName: Option[String] = dockerCredentials flatMap { _.authName }
  //val dockerToken: Option[String] = dockerCredentials map { _.token }

  val jobShell: String = configurationDescriptor.backendConfig.as[Option[String]]("job-shell").getOrElse(
    configurationDescriptor.globalConfig.getOrElse("system.job-shell", "/bin/bash"))
}
