package cromwell.backend.impl.jes

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.impl.jes.authentication.JesDockerCredentials
import cromwell.core.DockerConfiguration
import cromwell.filesystems.gcs.GoogleConfiguration

class JesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val root = configurationDescriptor.backendConfig.getString("root")
  val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)
  val jesAttributes = JesAttributes(googleConfig, configurationDescriptor.backendConfig)
  val dockerCredentials = DockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials map JesDockerCredentials.apply
  val needAuthFileUpload = jesAttributes.gcsFilesystemAuth.requiresAuthFile || dockerCredentials.isDefined
}
