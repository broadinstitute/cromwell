package cromwell.backend.google.pipelines.common

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.api.PipelinesApiFactoryInterface
import cromwell.backend.google.pipelines.common.authentication.PipelinesApiDockerCredentials
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.BackendDockerConfiguration


object PipelinesApiConfiguration {
  def apply(configurationDescriptor: BackendConfigurationDescriptor,
            genomicsFactory: PipelinesApiFactoryInterface,
            googleConfig: GoogleConfiguration,
            jesAttributes: PipelinesApiAttributes) = {
    new PipelinesApiConfiguration(configurationDescriptor, genomicsFactory, googleConfig, jesAttributes)
  }
}

class PipelinesApiConfiguration(
                                 val configurationDescriptor: BackendConfigurationDescriptor,
                                 val genomicsFactory: PipelinesApiFactoryInterface,
                                 val googleConfig: GoogleConfiguration,
                                 val jesAttributes: PipelinesApiAttributes,
                         ) {
  val jesAuths = jesAttributes.auths
  val root = configurationDescriptor.backendConfig.getString("root")
  val runtimeConfig = configurationDescriptor.backendRuntimeConfig
  val jesComputeServiceAccount = jesAttributes.computeServiceAccount
  val dockerCredentials = BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map PipelinesApiDockerCredentials.apply
  val needAuthFileUpload = jesAuths.gcs.requiresAuthFile || dockerCredentials.isDefined || jesAttributes.restrictMetadataAccess
  val qps = jesAttributes.qps
  val papiRequestWorkers = jesAttributes.requestWorkers
}
