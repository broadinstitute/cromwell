package cromwell.engine.backend.jes

import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.CreatePipelineRequest
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.Call
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.workflow.CallKey

import scala.collection.JavaConverters._

object Pipeline extends LazyLogging {
  def apply(command: String, workflow: WorkflowDescriptor, key: CallKey, jesParameters: Seq[JesParameter], projectId: String, jesConnection: JesInterface): Pipeline = {
    val call = key.scope
    val tag = s"JES Pipeline [UUID(${workflow.shortId}):${call.name}]"
    logger.debug(s"$tag Command line is: $command")
    val runtimeInfo = JesRuntimeInfo(command, call)

    val gcsPath = workflow.callDir(key)

    val cpr = new CreatePipelineRequest
    cpr.setProjectId(projectId)
    cpr.setDocker(runtimeInfo.docker)
    cpr.setResources(runtimeInfo.resources)
    cpr.setName(workflow.name)

    cpr.setParameters(jesParameters.map(_.toGoogleParameter).toVector.asJava)

    logger.info(s"$tag Pipeline parameters are:\n${cpr.getParameters.asScala.map(s=>s"  $s").mkString("\n")}")
    val pipelineId = jesConnection.genomics.pipelines().create(cpr).execute().getPipelineId
    logger.info(s"$tag Pipeline ID is $pipelineId")
    logger.info(s"$tag Project ID: $projectId")
    new Pipeline(command, 
                 pipelineId, 
                 projectId, 
                 gcsPath, 
                 workflow, 
                 call, 
                 jesParameters, 
                 runtimeInfo, 
                 jesConnection.genomics)
  }
}

// Note that id is the JES id not the workflow id
case class Pipeline(command: String,
                    id: String,
                    projectId: String,
                    gcsPath: String,
                    workflow: WorkflowDescriptor,
                    call: Call,
                    jesParameters: Seq[JesParameter],
                    runtimeInfo: JesRuntimeInfo,
                    genomicsService: Genomics) {
  def run: Run = Run(this)
}
