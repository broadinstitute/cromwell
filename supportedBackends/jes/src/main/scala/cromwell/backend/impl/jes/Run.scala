package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.labels.Labels
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
    val pipelineInfo = pipelineInfoBuilder.build(commandLine, runtimeAttributes)

    val pipeline = new Pipeline()
      .setProjectId(projectId)
      .setDocker(pipelineInfo.docker)
      .setResources(pipelineInfo.resources)
      .setName(workflow.workflow.unqualifiedName)
      .setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
      .setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    // disks cannot have mount points at runtime, so set them null
    val runtimePipelineResources = {
      val resources = pipelineInfoBuilder.build(commandLine, runtimeAttributes).resources
      val disksWithoutMountPoint = resources.getDisks.asScala map {
        _.setMountPoint(null)
      }
      resources.setDisks(disksWithoutMountPoint.asJava)
    }

    lazy val labels: Labels = {

      val subWorkflow = workflow.workflow
      val subWorkflowLabels = if (!subWorkflow.equals(workflow.rootWorkflow))
        Labels("cromwell-sub-workflow-name" -> subWorkflow.unqualifiedName)
      else
        Labels.empty

      val alias = jobDescriptor.call.unqualifiedName
      val aliasLabels = if (!alias.equals(jobDescriptor.call.task.name))
        Labels("wdl-call-alias" -> alias)
      else
        Labels.empty

      Labels(
        "cromwell-workflow-id" -> s"cromwell-${workflow.rootWorkflowId}",
        "cromwell-workflow-name" -> workflow.rootWorkflow.unqualifiedName,
        "wdl-task-name" -> jobDescriptor.call.task.name
      ) ++ subWorkflowLabels ++ aliasLabels ++ jobDescriptor.workflowDescriptor.customLabels
    }

    val svcAccount = new ServiceAccount().setEmail(computeServiceAccount).setScopes(GenomicsScopes)
    val rpargs = new RunPipelineArgs().setProjectId(projectId).setServiceAccount(svcAccount).setResources(runtimePipelineResources)

    rpargs.setInputs(jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
    rpargs.setOutputs(jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)

    rpargs.setLabels(labels.asJesLabels)

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
