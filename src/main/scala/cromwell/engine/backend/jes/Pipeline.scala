package cromwell.engine.backend.jes

import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.CreatePipelineRequest
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.workflow.CallKey
import cromwell.logging.WorkflowLogger
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._

object Pipeline {

  def apply(command: String,
            workflow: WorkflowDescriptor,
            key: CallKey,
            jesParameters: Seq[JesParameter],
            projectId: String,
            jesConnection: JesInterface,
            runIdForResumption: Option[String]): Pipeline = {

    val call = key.scope
    val logger = WorkflowLogger(
      "JES Pipeline",
      workflow,
      otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName)),
      callTag = Option(key.tag)
    )

    logger.debug(s"Command line is: $command")
    val runtimeInfo = JesRuntimeInfo(command, call)

    val gcsPath = workflow.callDir(key)

    val cpr = new CreatePipelineRequest
    cpr.setProjectId(projectId)
    cpr.setDocker(runtimeInfo.docker)
    cpr.setResources(runtimeInfo.resources)
    cpr.setName(workflow.name)

    cpr.setParameters(jesParameters.map(_.toGoogleParameter).toVector.asJava)

    def createPipeline = jesConnection.genomics.pipelines().create(cpr).execute().getPipelineId

    logger.info(s"Pipeline parameters are:\n${cpr.getParameters.asScala.map(s => s"  $s").mkString("\n")}")
    val pipelineId = if (runIdForResumption.isDefined) None else Option(createPipeline)

    logger.info(s"Pipeline ID is ${pipelineId.getOrElse("(none)")}")
    logger.info(s"Project ID: $projectId")
    new Pipeline(command,
                 pipelineId,
                 projectId,
                 gcsPath,
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
                            key: CallKey,
                            jesParameters: Seq[JesParameter],
                            runtimeInfo: JesRuntimeInfo,
                            genomicsService: Genomics,
                            runIdForResumption: Option[String]) {
  def run: Run = Run(this)
}
