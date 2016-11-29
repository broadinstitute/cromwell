package cromwell.backend.impl.jes

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.impl.jes.authentication.JesDockerCredentials
import cromwell.backend.impl.jes.io._
import cromwell.core.DockerConfiguration
import cromwell.core.path.CustomRetryParams
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.gcs.{GoogleConfiguration, RetryableGcsPathBuilderFactory}

import scala.concurrent.duration._
import scala.language.postfixOps

object JesConfiguration {
  val GcsRetryParams = CustomRetryParams(
    timeout = Duration.Inf,
    maxRetries = Option(3),
    backoff = SimpleExponentialBackoff(1 seconds, 3 seconds, 1.5D),
    isTransient = isTransientJesException,
    isFatal = isFatalJesException
  )
}

class JesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  private val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)

  val root = configurationDescriptor.backendConfig.getString("root")
  val jesAttributes = JesAttributes(googleConfig, configurationDescriptor.backendConfig)
  val jesAuths = jesAttributes.auths
  val jesComputeServiceAccount = jesAttributes.computeServiceAccount
  val gcsPathBuilderFactory = RetryableGcsPathBuilderFactory(jesAuths.gcs, customRetryParams = JesConfiguration.GcsRetryParams)
  val genomicsFactory = GenomicsFactory(googleConfig.applicationName, jesAuths.genomics, jesAttributes.endpointUrl)
  val dockerCredentials = DockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map JesDockerCredentials.apply
  val needAuthFileUpload = jesAuths.gcs.requiresAuthFile || dockerCredentials.isDefined
  val qps = jesAttributes.qps
}
