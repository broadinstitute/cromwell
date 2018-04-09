package cromwell.backend.impl.jes

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleFactory

// This is only a pointer to the implementation of Pipelines API v1 which is now in cromwell.backend.google.pipelines.v1alpha2
class JesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor) extends PipelinesApiLifecycleFactory(name, configurationDescriptor)
