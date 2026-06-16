package cromwell.backend.impl.tes

import akka.actor.Actor
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging

trait TesJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  lazy val initializationData: TesBackendInitializationData =
    backendInitializationDataAs[TesBackendInitializationData]

  lazy val tesWorkflowPaths: TesWorkflowPaths = workflowPaths.asInstanceOf[TesWorkflowPaths]

  lazy val tesJobPaths: TesJobPaths = jobPaths.asInstanceOf[TesJobPaths]

  lazy val tesConfiguration: TesConfiguration = initializationData.tesConfiguration

  lazy val runtimeAttributes = {
    val raw = jobDescriptor.runtimeAttributes
    jobLogger.debug(s"[backoff_limit] raw runtimeAttributes keys: ${raw.keys.mkString(", ")}")
    jobLogger.debug(s"[backoff_limit] raw runtimeAttributes values: ${raw.map { case (k, v) => s"$k=${v.toWomString}" }.mkString(", ")}")
    jobLogger.debug(s"[backoff_limit] useBackendParameters=${tesConfiguration.useBackendParameters}")
    val attrs = TesRuntimeAttributes(validatedRuntimeAttributes, raw, tesConfiguration)
    jobLogger.debug(s"[backoff_limit] resulting backendParameters: ${attrs.backendParameters}")
    attrs
  }
  override protected def nonStandardMetadata: Map[String, Any] =
    super.nonStandardMetadata ++ tesJobPaths.azureLogPathsForMetadata

}
