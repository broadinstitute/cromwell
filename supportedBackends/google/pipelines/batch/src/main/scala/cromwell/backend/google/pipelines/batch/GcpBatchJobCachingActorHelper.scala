package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging

trait GcpBatchJobCachingActorHelper extends StandardCachingActorHelper {
  this: GcpBatchAsyncBackendJobExecutionActor with JobLogging =>

  lazy val initializationData: GcpBackendInitializationData = {
    backendInitializationDataAs[GcpBackendInitializationData]
  }
  lazy val batchConfiguration: GcpBatchConfiguration = initializationData.gcpBatchConfiguration

  lazy val runtimeAttributes = GcpBatchRuntimeAttributes(
    validatedRuntimeAttributes,
    batchConfiguration
      .runtimeConfig
  )

  //lazy val workingDisk: PipelinesApiAttachedDisk = runtimeAttributes.disks.find(_.name == PipelinesApiWorkingDisk.Name).get

  //lazy val batchAttributes: GcpBatchConfigurationAttributes = batchConfiguration.papiAttributes

  //lazy val configuration: GcpBatchConfiguration = initializationData.configuration


}
