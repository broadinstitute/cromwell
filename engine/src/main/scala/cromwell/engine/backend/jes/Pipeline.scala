package cromwell.engine.backend.jes

import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model
import com.google.api.services.genomics.model.{Disk, PipelineParameter}
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.workflow.BackendCallKey
import cromwell.logging.WorkflowLogger
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._

object Pipeline {
  def apply(command: String,
            workflow: WorkflowDescriptor,
            key: BackendCallKey,
            runtimeAttributes: CromwellRuntimeAttributes,
            jesParameters: Seq[JesParameter],
            preemptible: Boolean,
            projectId: String,
            jesConnection: JesInterface,
            runIdForResumption: Option[String]): Pipeline = {
    val logger = WorkflowLogger(
      "JES Pipeline",
      workflow,
      otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName)),
      callTag = Option(key.tag)
    )

    logger.debug(s"Command line is: $command")
    val runtimeInfo = if (preemptible) PreemptibleJesRuntimeInfo(command, runtimeAttributes) else NonPreemptibleJesRuntimeInfo(command, runtimeAttributes)

    val gcsPath = workflow.callDir(key)

    val p = new model.Pipeline
    p.setProjectId(projectId)
    p.setDocker(runtimeInfo.docker)
    p.setResources(runtimeInfo.resources)
    p.setName(workflow.name)
    p.setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
    p.setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    def createPipeline = jesConnection.genomics.pipelines().create(p).execute().getPipelineId

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
    new Pipeline(command,
                 pipelineId,
                 projectId,
                 gcsPath.toString,
                 workflow,
                 key,
                 jesParameters,
                 runtimeInfo,
                 jesConnection.genomics,
                 runIdForResumption)
  }
}

case class Pipeline private(command: String,
                            pipelineId: Option[String],
                            projectId: String,
                            gcsPath: String,
                            workflow: WorkflowDescriptor,
                            key: BackendCallKey,
                            jesParameters: Seq[JesParameter],
                            runtimeInfo: JesRuntimeInfo,
                            genomicsService: Genomics,
                            runIdForResumption: Option[String]) {
  def run: Run = Run(this)
}
