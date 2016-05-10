package cromwell.backend.impl.jes

import com.google.api.services.genomics.model.{Disk, PipelineParameter}
import com.google.api.services.genomics.{Genomics, model}
import cromwell.backend.BackendJobDescriptor
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._

object Pipeline {
  def apply(jobDescriptor: BackendJobDescriptor,
            runtimeAttributes: JesRuntimeAttributes,
            callRootPath: String,
            commandLine: String,
            logFileName: String,
            jesParameters: Seq[JesParameter],
            projectId: String,
            genomicsInterface: Genomics,
            runIdForResumption: Option[String],
            preemptible: Boolean): Pipeline = {

    lazy val workflow = jobDescriptor.descriptor

    val logger = LoggerFactory.getLogger(Pipeline.getClass)

    logger.debug(s"Command line is: $commandLine")
    val runtimeInfoBuilder = if (preemptible) PreemptibleJesRuntimeInfoBuilder else NonPreemptibleJesRuntimeInfoBuilder
    val runtimeInfo = runtimeInfoBuilder.build(commandLine, runtimeAttributes)

    val p = new model.Pipeline
    p.setProjectId(projectId)
    p.setDocker(runtimeInfo.docker)
    p.setResources(runtimeInfo.resources)
    p.setName(workflow.workflowNamespace.workflow.unqualifiedName)
    p.setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
    p.setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    def createPipeline = genomicsInterface.pipelines().create(p).execute().getPipelineId

    def pipelineParameterString(p: PipelineParameter): String = {
      val description = Option(p.getLocalCopy) match {
        case Some(localCopy) => s"disk:${localCopy.getDisk} relpath:${localCopy.getPath}"
        case None => s"(literal value)"
      }
      s"  ${p.getName} -> $description"
    }
    def diskString(d: Disk): String = {
      s"  ${d.getName} -> ${d.getMountPoint} (${d.getSizeGb}GB ${d.getType})"
    }
    logger.info(s"Inputs:\n${p.getInputParameters.asScala.map(pipelineParameterString).mkString("\n")}")
    logger.info(s"Outputs:\n${p.getOutputParameters.asScala.map(pipelineParameterString).mkString("\n")}")
    logger.info(s"Mounts:\n${runtimeInfo.resources.getDisks.asScala.map(diskString).mkString("\n")}")

    val pipelineId = if (runIdForResumption.isDefined) None else Option(createPipeline)

    logger.info(s"Pipeline ID is ${pipelineId.getOrElse("(none)")}")
    logger.info(s"Project ID: $projectId")
    new Pipeline(
      jobDescriptor = jobDescriptor,
      command = commandLine,
      logFileName = logFileName,
      pipelineId = pipelineId,
      projectId = projectId,
      gcsPath = callRootPath,
      jesParameters = jesParameters,
      runtimeInfo = runtimeInfo,
      genomicsService = genomicsInterface)
  }
}

case class Pipeline private(jobDescriptor: BackendJobDescriptor,
                            logFileName: String,
                            command: String,
                            pipelineId: Option[String],
                            projectId: String,
                            gcsPath: String,
                            jesParameters: Seq[JesParameter],
                            runtimeInfo: JesRuntimeInfo,
                            genomicsService: Genomics,
                            runIdForResumption: Option[String] = None) {
  def run: Run = Run(this, logFileName)
}
