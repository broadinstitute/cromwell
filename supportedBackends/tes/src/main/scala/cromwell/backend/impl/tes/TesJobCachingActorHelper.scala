package cromwell.backend.impl.tes


import akka.actor.Actor
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging

trait TesJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  lazy val initializationData: TesBackendInitializationData = {
    backendInitializationDataAs[TesBackendInitializationData]
  }

  lazy val tesWorkflowPaths: TesWorkflowPaths = workflowPaths.asInstanceOf[TesWorkflowPaths]

  lazy val tesJobPaths: TesJobPaths = jobPaths.asInstanceOf[TesJobPaths]

  lazy val tesConfiguration: TesConfiguration = initializationData.tesConfiguration

  lazy val runtimeAttributes = TesRuntimeAttributes(validatedRuntimeAttributes, tesConfiguration.runtimeConfig)
}
