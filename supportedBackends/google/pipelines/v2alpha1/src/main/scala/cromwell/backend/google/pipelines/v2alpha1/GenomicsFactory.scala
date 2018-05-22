package cromwell.backend.google.pipelines.v2alpha1

import java.net.URL

import com.google.api.client.http.HttpRequestInitializer
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.{ActionBuilder, Delocalization, Localization}
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.StandardAsyncJob
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

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

    override def runRequest(createPipelineParameters: CreatePipelineParameters, jobPaths: Option[JobPaths]) = {
      // Disks defined in the runtime attributes
      val disks = createPipelineParameters.toDisks
      // Mounts for disks defined in the runtime attributes
      val mounts = createPipelineParameters.toMounts

      val localization: List[Action] = localizeActions(createPipelineParameters, mounts)
      // localization.size + 1 because action indices are 1-based and the next action after localization will be the user's
      val deLocalization: List[Action] = deLocalizeActions(createPipelineParameters, mounts, localization.size + 1, jobPaths)

      val environment = Map.empty[String, String].asJava

      val userAction = ActionBuilder.userAction(
        createPipelineParameters.dockerImage,
        createPipelineParameters.commandScriptContainerPath.pathAsString,
        mounts
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
        .gpuResource.map(_.toAccelerator).toList.asJava

      val virtualMachine = new VirtualMachine()
        .setDisks(disks.asJava)
        .setPreemptible(createPipelineParameters.preemptible)
        .setServiceAccount(serviceAccount)
        .setMachineType(createPipelineParameters.runtimeAttributes.toMachineType)
        .setBootDiskSizeGb(createPipelineParameters.runtimeAttributes.bootDiskSize)
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
