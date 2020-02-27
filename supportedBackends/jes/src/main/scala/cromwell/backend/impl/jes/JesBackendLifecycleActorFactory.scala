package cromwell.backend.impl.jes

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory
import org.slf4j.LoggerFactory

// This is only a pointer to the implementation of Pipelines API v1 which is now in cromwell.backend.google.pipelines.v1alpha2
class JesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor) extends PipelinesApiLifecycleActorFactory(name, configurationDescriptor) {
  val logger = LoggerFactory.getLogger("DeprecatedJesBackend")
  logger.warn("This actor factory is deprecated. Please use cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory for PAPI v1 or cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory for PAPI v2")
}
