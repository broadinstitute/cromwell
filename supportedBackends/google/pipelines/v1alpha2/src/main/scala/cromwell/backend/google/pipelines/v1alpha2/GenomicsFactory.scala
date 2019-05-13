package cromwell.backend.google.pipelines.v1alpha2

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.genomics.{Genomics, GenomicsScopes}
import com.google.api.services.genomics.model._
import common.collections.EnhancedCollections._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.common.{PipelinesApiInput, PipelinesApiOutput}
import cromwell.backend.google.pipelines.v1alpha2.PipelinesConversions._
import cromwell.backend.google.pipelines.v1alpha2.ToParameter.ops._
import cromwell.backend.google.pipelines.v1alpha2.api.{NonPreemptiblePipelineInfoBuilder, PreemptiblePipelineInfoBuilder}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.logging.JobLogger

import scala.collection.JavaConverters._

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) extends PipelinesApiFactoryInterface {
  def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = {
    new PipelinesApiRequestFactory {
      private val genomics = new Genomics.Builder(
        GoogleAuthMode.httpTransport,
        GoogleAuthMode.jsonFactory,
        initializer)
        .setApplicationName(applicationName)
        .setRootUrl(endpointUrl.toString)
        .build

      override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger) = {
        lazy val workflow = createPipelineParameters.jobDescriptor.workflowDescriptor
        val commandLine = s"/bin/bash ${createPipelineParameters.commandScriptContainerPath.pathAsString}"
        
        val pipelineInfoBuilder = if (createPipelineParameters.preemptible) PreemptiblePipelineInfoBuilder else NonPreemptiblePipelineInfoBuilder
        val pipelineInfo = pipelineInfoBuilder.build(commandLine, createPipelineParameters.runtimeAttributes, createPipelineParameters.dockerImage)

        val inputParameters: Map[String, PipelinesApiInput] = createPipelineParameters.inputParameters.map({ i => i.name -> i }).toMap
        val outputParameters: Map[String, PipelinesApiOutput] = createPipelineParameters.outputParameters.map({ o => o.name -> o }).toMap

        val literalInputPipelineParamters = createPipelineParameters.literalInputs.map(_.toGooglePipelineParameter)
        val literalInputRunParameters = createPipelineParameters.literalInputs.map(l => l.name -> l.toGoogleRunParameter).toMap

        val pipeline = new Pipeline()
          .setProjectId(createPipelineParameters.projectId)
          .setDocker(pipelineInfo.docker)
          .setResources(pipelineInfo.resources)
          .setName(workflow.callable.name)
          .setInputParameters((inputParameters.values.map(_.toGooglePipelineParameter) ++ literalInputPipelineParamters).toVector.asJava)
          .setOutputParameters(outputParameters.values.map(_.toGooglePipelineParameter).toVector.asJava)

        // disks cannot have mount points at runtime, so set them null
        val runtimePipelineResources = {
          val resources = pipelineInfoBuilder.build(commandLine, createPipelineParameters.runtimeAttributes, createPipelineParameters.dockerImage).resources
          val disksWithoutMountPoint = resources.getDisks.asScala map {
            _.setMountPoint(null)
          }
          resources.setDisks(disksWithoutMountPoint.asJava)
        }

        val svcAccount = new ServiceAccount().setEmail(createPipelineParameters.computeServiceAccount)
          .setScopes(List(
            GenomicsScopes.GENOMICS,
            ComputeScopes.COMPUTE
          ).asJava)
        val rpargs = new RunPipelineArgs().setProjectId(createPipelineParameters.projectId).setServiceAccount(svcAccount).setResources(runtimePipelineResources)

        rpargs.setInputs((inputParameters.safeMapValues(_.toGoogleRunParameter) ++ literalInputRunParameters).asJava)
        rpargs.setOutputs(outputParameters.safeMapValues(_.toGoogleRunParameter).asJava)

        rpargs.setLabels(createPipelineParameters.googleLabels.map(label => label.key -> label.value).toMap.asJava)

        val rpr = new RunPipelineRequest().setEphemeralPipeline(pipeline).setPipelineArgs(rpargs)

        val logging = new LoggingOptions()
        logging.setGcsPath(createPipelineParameters.logGcsPath.pathAsString)
        rpargs.setLogging(logging)

        genomics.pipelines().run(rpr).buildHttpRequest()
      }

      override def getRequest(job: StandardAsyncJob) = genomics.operations().get(job.jobId).buildHttpRequest()

      override def cancelRequest(job: StandardAsyncJob): HttpRequest = {
        val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
        genomics.operations().cancel(job.jobId, cancellationRequest).buildHttpRequest()
      }
    }
  }

  override def usesEncryptedDocker: Boolean = false
}
