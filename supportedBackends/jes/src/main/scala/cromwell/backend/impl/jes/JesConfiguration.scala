package cromwell.backend.impl.jes

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.impl.jes.authentication.JesDockerCredentials
import cromwell.core.BackendDockerConfiguration
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}

class JesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)

  val root = configurationDescriptor.backendConfig.getString("root")
  val runtimeConfig = configurationDescriptor.backendRuntimeConfig
  val jesAttributes = JesAttributes(googleConfig, configurationDescriptor.backendConfig)
  val jesAuths = jesAttributes.auths
  val jesComputeServiceAccount = jesAttributes.computeServiceAccount
  val gcsPathBuilderFactory = GcsPathBuilderFactory(jesAuths.gcs, googleConfig.applicationName)
  val genomicsFactory = GenomicsFactory(googleConfig.applicationName, jesAuths.genomics, jesAttributes.endpointUrl)
  val dockerCredentials = BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map JesDockerCredentials.apply
  val needAuthFileUpload = jesAuths.gcs.requiresAuthFile || dockerCredentials.isDefined
  val qps = jesAttributes.qps
}
