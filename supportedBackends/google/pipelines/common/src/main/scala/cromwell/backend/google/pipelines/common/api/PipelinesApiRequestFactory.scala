package cromwell.backend.google.pipelines.common.api

import com.google.api.client.http.HttpRequest
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.{JesParameter, PipelinesApiRuntimeAttributes}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.labels.Labels

trait PipelinesApiRequestFactory {
  // Genomics requests
  def makeRunPipelineRequest(createPipelineParameters: CreatePipelineParameters): HttpRequest
  def getOperationRequest(job: StandardAsyncJob): HttpRequest
  def abortRequest(job: StandardAsyncJob): HttpRequest
}

object PipelinesApiRequestFactory {
  case class CreatePipelineParameters(jobDescriptor: BackendJobDescriptor,
                                      runtimeAttributes: PipelinesApiRuntimeAttributes,
                                      dockerImage: String,
                                      callRootPath: String,
                                      commandLine: String,
                                      logFileName: String,
                                      jesParameters: Seq[JesParameter],
                                      projectId: String,
                                      computeServiceAccount: String,
                                      labels: Labels,
                                      preemptible: Boolean)
}
