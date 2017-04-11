package cromwell.backend.impl.jes

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.impl.jes.authentication.JesDockerCredentials
import cromwell.core.DockerConfiguration
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}

class JesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)

  val root = configurationDescriptor.backendConfig.getString("root")
  val jesAttributes = JesAttributes(googleConfig, configurationDescriptor.backendConfig)
  val jesAuths = jesAttributes.auths
  val jesComputeServiceAccount = jesAttributes.computeServiceAccount
  val gcsPathBuilderFactory = GcsPathBuilderFactory(jesAuths.gcs, googleConfig.applicationName)
  val genomicsFactory = GenomicsFactory(googleConfig.applicationName, jesAuths.genomics, jesAttributes.endpointUrl)
  val dockerCredentials = DockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map JesDockerCredentials.apply
  val needAuthFileUpload = jesAuths.gcs.requiresAuthFile || dockerCredentials.isDefined
  val qps = jesAttributes.qps
  val defaultZones = jesAttributes.defaultZones
}
