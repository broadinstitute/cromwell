package cromwell.backend.google.batch.api

import com.google.cloud.batch.v1.{CreateJobRequest, DeleteJobRequest, GetJobRequest, JobName}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.batch.io.GcpBatchAttachedDisk
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.VirtualPrivateCloudConfiguration
import cromwell.backend.google.batch.models._
import cromwell.backend.google.batch.monitoring.{CheckpointingConfiguration, MonitoringImage}
import cromwell.core.logging.JobLogger
import cromwell.core.path.Path
import wom.runtime.WomOutputRuntimeExtractor

import scala.concurrent.duration.FiniteDuration

trait GcpBatchRequestFactory {
  def submitRequest(data: GcpBatchRequest, jobLogger: JobLogger): CreateJobRequest

  def queryRequest(jobName: JobName): GetJobRequest

  def abortRequest(jobName: JobName): DeleteJobRequest

}

object GcpBatchRequestFactory {

  type MountsToEnv = List[String] => Map[String, String]

  /**
    * Input parameters that are not strictly needed by the user's command but are Cromwell byproducts.
    */
  case class DetritusInputParameters(
    executionScriptInputParameter: GcpBatchFileInput,
    monitoringScriptInputParameter: Option[GcpBatchFileInput]
  ) {
    def all: List[GcpBatchFileInput] = List(executionScriptInputParameter) ++ monitoringScriptInputParameter
  }

  /**
    * Output parameters that are not produced by the user's command but are Cromwell byproducts.
    */
  case class DetritusOutputParameters(
    monitoringScriptOutputParameter: Option[GcpBatchFileOutput],
    rcFileOutputParameter: GcpBatchFileOutput,
    memoryRetryRCFileOutputParameter: GcpBatchFileOutput,
    logFileOutputParameter: GcpBatchFileOutput
  ) {
    def all: List[GcpBatchFileOutput] = memoryRetryRCFileOutputParameter ::
      logFileOutputParameter ::
      rcFileOutputParameter ::
      monitoringScriptOutputParameter.toList
  }

  /**
    * Bundle containing all input and output parameters to a Batch job
    * Detrituses and actual inputs / outputs to the job are separated for more clarity and to leave open the possibility
    * to treat them differently.
    */
  case class InputOutputParameters(
    detritusInputParameters: DetritusInputParameters,
    jobInputParameters: List[GcpBatchInput],
    jobOutputParameters: List[GcpBatchOutput],
    detritusOutputParameters: DetritusOutputParameters
  ) {
    lazy val fileInputParameters: List[GcpBatchInput] = jobInputParameters ++ detritusInputParameters.all
    lazy val fileOutputParameters: List[GcpBatchOutput] = detritusOutputParameters.all ++ jobOutputParameters
  }

  case class CreateBatchDockerKeyAndToken(key: String, encryptedToken: String)

  case class CreateBatchJobParameters(jobDescriptor: BackendJobDescriptor,
                                      runtimeAttributes: GcpBatchRuntimeAttributes,
                                      dockerImage: String,
                                      cloudWorkflowRoot: Path,
                                      cloudCallRoot: Path,
                                      commandScriptContainerPath: Path,
                                      logGcsPath: Path,
                                      inputOutputParameters: InputOutputParameters,
                                      projectId: String,
                                      computeServiceAccount: String,
                                      googleLabels: Seq[GcpLabel],
                                      batchTimeout: FiniteDuration,
                                      jobShell: String,
                                      privateDockerKeyAndEncryptedToken: Option[CreateBatchDockerKeyAndToken],
                                      womOutputRuntimeExtractor: Option[WomOutputRuntimeExtractor],
                                      disks: Seq[GcpBatchAttachedDisk],
                                      virtualPrivateCloudConfiguration: VirtualPrivateCloudConfiguration,
                                      retryWithMoreMemoryKeys: Option[List[String]],
                                      fuseEnabled: Boolean,
                                      referenceDisksForLocalizationOpt: Option[List[GcpBatchAttachedDisk]],
                                      monitoringImage: MonitoringImage,
                                      checkpointingConfiguration: CheckpointingConfiguration,
                                      enableSshAccess: Boolean,
                                      vpcNetworkAndSubnetworkProjectLabels: Option[VpcAndSubnetworkProjectLabelValues],
                                      dockerhubCredentials: (String, String)
  ) {
    def outputParameters = inputOutputParameters.fileOutputParameters
  }

}
