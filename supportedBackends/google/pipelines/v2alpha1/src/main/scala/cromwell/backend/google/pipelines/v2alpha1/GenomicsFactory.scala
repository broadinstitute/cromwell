package cromwell.backend.google.pipelines.v2alpha1

import java.net.URL

import com.google.api.client.http.HttpRequestInitializer
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.{ActionBuilder, Delocalization, Localization}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.logging.JobLogger
import mouse.all._

import scala.collection.JavaConverters._

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) extends PipelinesApiFactoryInterface
  with Localization
  with Delocalization {

  override def build(initializer: HttpRequestInitializer): PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
    val genomics = new Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      initializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def cancelRequest(job: StandardAsyncJob) = {
      genomics.projects().operations().cancel(job.jobId, new CancelOperationRequest()).buildHttpRequest()
    }

    override def getRequest(job: StandardAsyncJob) = {
      genomics.projects().operations().get(job.jobId).buildHttpRequest()
    }

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobLogger: JobLogger) = {
      // Disks defined in the runtime attributes
      val disks = createPipelineParameters |> toDisks
      // Mounts for disks defined in the runtime attributes
      val mounts = createPipelineParameters |> toMounts

      val localization: List[Action] = localizeActions(createPipelineParameters, mounts)
      val deLocalization: List[Action] = deLocalizeActions(createPipelineParameters, mounts)

      val environment = Map.empty[String, String].asJava

      val userAction = ActionBuilder.userAction(
        createPipelineParameters.dockerImage,
        createPipelineParameters.commandScriptContainerPath.pathAsString,
        mounts,
        createPipelineParameters.jobShell
      )

      val serviceAccount = new ServiceAccount()
        .setEmail(createPipelineParameters.computeServiceAccount)
        .setScopes(
          List(
            PipelinesApiFactoryInterface.GenomicsScope,
            PipelinesApiFactoryInterface.ComputeScope,
            PipelinesApiFactoryInterface.StorageFullControlScope
          ).asJava
        )

      val network = new Network()
        .setUsePrivateAddress(createPipelineParameters.runtimeAttributes.noAddress)

      val accelerators = createPipelineParameters.runtimeAttributes
        .gpuResource.map(toAccelerator).toList.asJava
      
      /*
       * Adjust for the additional docker images Cromwell uses for (de)localization
       * At the moment, google/cloud-sdk:alpine (173MB) and stedolan/jq:latest (182MB)
       * Round it up to 1GB
       */
      val adjustedBootDiskSize = createPipelineParameters.runtimeAttributes.bootDiskSize + 1

      val virtualMachine = new VirtualMachine()
        .setDisks(disks.asJava)
        .setPreemptible(createPipelineParameters.preemptible)
        .setServiceAccount(serviceAccount)
        .setMachineType(createPipelineParameters.runtimeAttributes |> toMachineType(jobLogger))
        .setBootDiskSizeGb(adjustedBootDiskSize)
        .setLabels(createPipelineParameters.labels.asJavaMap)
        .setNetwork(network)
        .setAccelerators(accelerators)

      createPipelineParameters.runtimeAttributes.gpuResource foreach { resource =>
        virtualMachine.setNvidiaDriverVersion(resource.nvidiaDriverVersion)
      }

      val resources = new Resources()
        .setProjectId(createPipelineParameters.projectId)
        .setZones(createPipelineParameters.runtimeAttributes.zones.asJava)
        .setVirtualMachine(virtualMachine)

      val pipeline = new Pipeline()
        .setResources(resources)
        .setActions((localization ++ List(userAction) ++ deLocalization).asJava)
        .setEnvironment(environment)

      val pipelineRequest = new RunPipelineRequest()
        .setPipeline(pipeline)
        .setLabels(createPipelineParameters.labels.asJavaMap)

      genomics.pipelines().run(pipelineRequest).buildHttpRequest()
    }
  }
}
