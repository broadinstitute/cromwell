package cromwell.backend.google.pipelines.v1alpha2

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.v1alpha2.GenomicsFactory._
import cromwell.backend.google.pipelines.v1alpha2.PipelinesImplicits._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import scala.collection.JavaConverters._

object GenomicsFactory {
  private val GenomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava
}

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) extends PipelinesApiFactoryInterface {
  def fromCredentials(credentials: Credentials): PipelinesApiRequestFactory = {
    val httpRequestInitializer = {
      val delegate = new HttpCredentialsAdapter(credentials)
      new HttpRequestInitializer() {
        def initialize(httpRequest: HttpRequest) = {
          delegate.initialize(httpRequest)
        }
      }
    }

    new PipelinesApiRequestFactory {
      private val genomics = new Genomics.Builder(
        GoogleAuthMode.httpTransport,
        GoogleAuthMode.jsonFactory,
        httpRequestInitializer)
        .setApplicationName(applicationName)
        .setRootUrl(endpointUrl.toString)
        .build

      override def makeRunPipelineRequest(createPipelineParameters: CreatePipelineParameters) = {
        lazy val workflow = createPipelineParameters.jobDescriptor.workflowDescriptor
        val pipelineInfoBuilder = if (createPipelineParameters.preemptible) PreemptibleJesPipelineInfoBuilder else NonPreemptibleJesPipelineInfoBuilder
        val pipelineInfo = pipelineInfoBuilder.build(createPipelineParameters.commandLine, createPipelineParameters.runtimeAttributes, createPipelineParameters.dockerImage)

        val pipeline = new Pipeline()
          .setProjectId(createPipelineParameters.projectId)
          .setDocker(pipelineInfo.docker)
          .setResources(pipelineInfo.resources)
          .setName(workflow.callable.name)
          .setInputParameters(createPipelineParameters.jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
          .setOutputParameters(createPipelineParameters.jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

        // disks cannot have mount points at runtime, so set them null
        val runtimePipelineResources = {
          val resources = pipelineInfoBuilder.build(createPipelineParameters.commandLine, createPipelineParameters.runtimeAttributes, createPipelineParameters.dockerImage).resources
          val disksWithoutMountPoint = resources.getDisks.asScala map {
            _.setMountPoint(null)
          }
          resources.setDisks(disksWithoutMountPoint.asJava)
        }

        val svcAccount = new ServiceAccount().setEmail(createPipelineParameters.computeServiceAccount).setScopes(GenomicsScopes)
        val rpargs = new RunPipelineArgs().setProjectId(createPipelineParameters.projectId).setServiceAccount(svcAccount).setResources(runtimePipelineResources)

        rpargs.setInputs(createPipelineParameters.jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
        rpargs.setOutputs(createPipelineParameters.jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)

        rpargs.setLabels(createPipelineParameters.labels.asJavaMap)

        val rpr = new RunPipelineRequest().setEphemeralPipeline(pipeline).setPipelineArgs(rpargs)

        val logging = new LoggingOptions()
        logging.setGcsPath(s"${createPipelineParameters.callRootPath}/${createPipelineParameters.logFileName}")
        rpargs.setLogging(logging)

        genomics.pipelines().run(rpr).buildHttpRequest()
      }

      override def getOperationRequest(job: StandardAsyncJob) = genomics.operations().get(job.jobId).buildHttpRequest()

      override def abortRequest(job: StandardAsyncJob): HttpRequest = {
        val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
        genomics.operations().cancel(job.jobId, cancellationRequest).buildHttpRequest()
      }
    }
  }
}
