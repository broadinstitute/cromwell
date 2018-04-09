package cromwell.backend.google.pipelines.v2alpha1

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.{PipelinesApiBackendLifecycleActorFactory, PipelinesApiConfiguration}

class PipelinesApiLifecycleFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends PipelinesApiBackendLifecycleActorFactory(name, configurationDescriptor) {
  override protected val jesConfiguration = PipelinesApiConfiguration(configurationDescriptor, null, googleConfig, jesAttributes)
}
