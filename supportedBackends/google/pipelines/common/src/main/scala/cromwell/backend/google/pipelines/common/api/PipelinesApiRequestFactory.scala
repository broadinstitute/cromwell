package cromwell.backend.google.pipelines.common.api

import com.google.api.client.http.HttpRequest
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.labels.Labels
import cromwell.core.path.Path

/**
  * The PipelinesApiRequestFactory defines the HttpRequests needed to run jobs
  */
trait PipelinesApiRequestFactory {
  def runRequest(createPipelineParameters: CreatePipelineParameters, jobPaths: Option[JobPaths]): HttpRequest
  def getRequest(job: StandardAsyncJob): HttpRequest
  def cancelRequest(job: StandardAsyncJob): HttpRequest
}

object PipelinesApiRequestFactory {

  /**
    * Input parameters that are not strictly needed by the user's command but are Cromwell byproducts.
    */
  case class DetritusInputParameters(
                                      executionScriptInputParameter: PipelinesApiFileInput,
                                      monitoringScriptInputParameter: Option[PipelinesApiFileInput]
                                    ) {
    def all: List[PipelinesApiFileInput] = List(executionScriptInputParameter) ++ monitoringScriptInputParameter
  }

  /**
    * Output parameters that are not produced by the user's command but are Cromwell byproducts.
    */
  case class DetritusOutputParameters(
                                       monitoringScriptOutputParameter: Option[PipelinesApiFileOutput],
                                       rcFileOutputParameter: PipelinesApiFileOutput
                                    ) {
    def all: List[PipelinesApiFileOutput] = List(rcFileOutputParameter) ++ monitoringScriptOutputParameter
  }

  /**
    * Bundle containing all input and output parameters to a PAPI job
    * Detrituses and actual inputs / outputs to the job are separated for more clarity and to leave open the possibility
    * to treat them differently.
    */
  case class InputOutputParameters(
                                    detritusInputParameters: DetritusInputParameters,
                                    jobInputParameters: List[PipelinesApiFileInput],
                                    jobOutputParameters: List[PipelinesApiFileOutput],
                                    detritusOutputParameters: DetritusOutputParameters,
                                    literalInputParameters: List[PipelinesApiLiteralInput]
                                  ) {
    lazy val fileInputParameters: List[PipelinesApiFileInput] = jobInputParameters ++ detritusInputParameters.all
    lazy val fileOutputParameters: List[PipelinesApiFileOutput] = jobOutputParameters ++ detritusOutputParameters.all
  }
  
  case class CreatePipelineParameters(jobDescriptor: BackendJobDescriptor,
                                      runtimeAttributes: PipelinesApiRuntimeAttributes,
                                      dockerImage: String,
                                      cloudCallRoot: Path,
                                      commandScriptContainerPath: Path,
                                      logGcsPath: Path,
                                      inputOutputParameters: InputOutputParameters,
                                      projectId: String,
                                      computeServiceAccount: String,
                                      labels: Labels,
                                      preemptible: Boolean) {
    def inputParameters = inputOutputParameters.literalInputParameters ++ inputOutputParameters.fileInputParameters
    def outputParameters = inputOutputParameters.fileOutputParameters
    def allParameters = inputParameters ++ outputParameters
  }
}
