package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path

trait GcpBatchJobCachingActorHelper extends StandardCachingActorHelper {
  this: GcpBatchAsyncBackendJobExecutionActor with JobLogging =>

  lazy val initializationData: GcpBackendInitializationData = {
    backendInitializationDataAs[GcpBackendInitializationData]
  }
  lazy val batchConfiguration: GcpBatchConfiguration = initializationData.gcpBatchConfiguration

  lazy val gcpBatchCallPaths: GcpBatchJobPaths = jobPaths.asInstanceOf[GcpBatchJobPaths]

  lazy val runtimeAttributes = GcpBatchRuntimeAttributes(
    validatedRuntimeAttributes,
    batchConfiguration
      .runtimeConfig
  )

  //lazy val workingDisk: PipelinesApiAttachedDisk = runtimeAttributes.disks.find(_.name == PipelinesApiWorkingDisk.Name).get

  lazy val callRootPath: Path = gcpBatchCallPaths.callExecutionRoot
  lazy val returnCodeFilename: String = gcpBatchCallPaths.returnCodeFilename
  lazy val returnCodeGcsPath: Path = gcpBatchCallPaths.returnCode
  //lazy val jesLogPath: Path = gcpBatchCallPaths.
  lazy val memoryRetryRCFilename: String = gcpBatchCallPaths.memoryRetryRCFilename
  lazy val memoryRetryRCGcsPath: Path = gcpBatchCallPaths.memoryRetryRC

  lazy val batchAttributes: GcpBatchConfigurationAttributes = batchConfiguration.batchAttributes

  //lazy val configuration: GcpBatchConfiguration = initializationData.configuration


}
