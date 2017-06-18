package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.standard.StandardAsyncJob
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._

object Run {
  private val GenomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava

  val NoAddressFieldName = "noAddress"

  val slf4jLogger = LoggerFactory.getLogger(Run.getClass)

  def makeRunPipelineRequest(jobDescriptor: BackendJobDescriptor,
                             runtimeAttributes: JesRuntimeAttributes,
                             dockerImage: String,
                             callRootPath: String,
                             commandLine: String,
                             logFileName: String,
                             jesParameters: Seq[JesParameter],
                             projectId: String,
                             computeServiceAccount: String,
                             preemptible: Boolean,
                             genomicsInterface: Genomics): RunPipelineRequest = {

    lazy val workflow = jobDescriptor.workflowDescriptor
    val pipelineInfoBuilder = if (preemptible) PreemptibleJesPipelineInfoBuilder else NonPreemptibleJesPipelineInfoBuilder
    val pipelineInfo = pipelineInfoBuilder.build(commandLine, runtimeAttributes, dockerImage)

    val pipeline = new Pipeline()
      .setProjectId(projectId)
      .setDocker(pipelineInfo.docker)
      .setResources(pipelineInfo.resources)
      .setName(workflow.workflow.unqualifiedName)
      .setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
      .setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    // disks cannot have mount points at runtime, so set them null
    val runtimePipelineResources = {
      val resources = pipelineInfoBuilder.build(commandLine, runtimeAttributes, dockerImage).resources
      val disksWithoutMountPoint = resources.getDisks.asScala map {
        _.setMountPoint(null)
      }
      resources.setDisks(disksWithoutMountPoint.asJava)
    }

    val svcAccount = new ServiceAccount().setEmail(computeServiceAccount).setScopes(GenomicsScopes)
    val rpargs = new RunPipelineArgs().setProjectId(projectId).setServiceAccount(svcAccount).setResources(runtimePipelineResources)

    rpargs.setInputs(jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
    rpargs.setOutputs(jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)

    rpargs.setLabels(workflow.customLabels.asJavaMap)

    val rpr = new RunPipelineRequest().setEphemeralPipeline(pipeline).setPipelineArgs(rpargs)

    val logging = new LoggingOptions()
    logging.setGcsPath(s"$callRootPath/$logFileName")
    rpargs.setLogging(logging)

    rpr
  }
}

final case class Run(job: StandardAsyncJob, genomicsInterface: Genomics) {

  def getOperationCommand = genomicsInterface.operations().get(job.jobId)

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    genomicsInterface.operations().cancel(job.jobId, cancellationRequest).execute
    ()
  }
}
