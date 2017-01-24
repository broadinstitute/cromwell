package cromwell.backend.impl.tes


import akka.actor.Actor
import cromwell.backend.io.{JobPathsWithDocker, WorkflowPathsWithDocker}
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging

trait TesJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  lazy val initializationData: TesBackendInitializationData = {
    backendInitializationDataAs[TesBackendInitializationData]
  }

  lazy val tesWorkflowPaths: WorkflowPathsWithDocker = workflowPaths.asInstanceOf[WorkflowPathsWithDocker]

  lazy val tesJobPaths: JobPathsWithDocker = jobPaths.asInstanceOf[JobPathsWithDocker]

  lazy val tesConfiguration: TesConfiguration = initializationData.tesConfiguration

  lazy val runtimeAttributes = TesRuntimeAttributes(validatedRuntimeAttributes)
}
