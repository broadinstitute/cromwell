package cromwell.backend.google.pipelines.v1alpha2

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import cromwell.backend.google.pipelines.common.{PipelinesApiFileOutput, PipelinesApiInput}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.v1alpha2.PipelinesConversions._
import cromwell.backend.google.pipelines.v1alpha2.ToParameter.ops._
import cromwell.backend.google.pipelines.v1alpha2.api.{NonPreemptiblePipelineInfoBuilder, PreemptiblePipelineInfoBuilder}
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

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

      override def runRequest(createPipelineParameters: CreatePipelineParameters, jobPaths: Option[JobPaths]) = {
        lazy val workflow = createPipelineParameters.jobDescriptor.workflowDescriptor
        val commandLine = s"/bin/bash ${createPipelineParameters.commandScriptContainerPath.pathAsString}"
        
        val pipelineInfoBuilder = if (createPipelineParameters.preemptible) PreemptiblePipelineInfoBuilder else NonPreemptiblePipelineInfoBuilder
        val pipelineInfo = pipelineInfoBuilder.build(commandLine, createPipelineParameters.runtimeAttributes, createPipelineParameters.dockerImage)

        val inputParameters: Map[String, PipelinesApiInput] = createPipelineParameters.inputParameters.map({ i => i.name -> i }).toMap
        val outputParameters: Map[String, PipelinesApiFileOutput] = createPipelineParameters.outputParameters.map({ o => o.name -> o }).toMap

        val pipeline = new Pipeline()
          .setProjectId(createPipelineParameters.projectId)
          .setDocker(pipelineInfo.docker)
          .setResources(pipelineInfo.resources)
          .setName(workflow.callable.name)
          .setInputParameters(inputParameters.values.map(_.toGooglePipelineParameter).toVector.asJava)
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
            PipelinesApiFactoryInterface.GenomicsScope,
            PipelinesApiFactoryInterface.ComputeScope
          ).asJava)
        val rpargs = new RunPipelineArgs().setProjectId(createPipelineParameters.projectId).setServiceAccount(svcAccount).setResources(runtimePipelineResources)

        rpargs.setInputs(inputParameters.mapValues(_.toGoogleRunParameter).asJava)
        rpargs.setOutputs(outputParameters.mapValues(_.toGoogleRunParameter).asJava)

        rpargs.setLabels(createPipelineParameters.labels.asJavaMap)

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
}
