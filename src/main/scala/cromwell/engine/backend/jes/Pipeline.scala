package cromwell.engine.backend.jes


import cromwell.binding.{Call, WorkflowDescriptor}
import scala.collection.JavaConverters._
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.CreatePipelineRequest
import JesBackend._



object Pipeline {
  def apply(command: String, workflow: WorkflowDescriptor, call: Call, jesParameters: Seq[JesParameter], projectId: String, jesConnection: JesInterface): Pipeline = {
    println(s"Command line is $command")
    val runtimeInfo = JesRuntimeInfo(command, call)

    val gcsPath = workflow.callDir(call)

    val cpr = new CreatePipelineRequest
    cpr.setProjectId(projectId)
    cpr.setDocker(runtimeInfo.docker)
    cpr.setResources(runtimeInfo.resources) // FIXME: These are still F'd up if you dig into it
    cpr.setName(workflow.name)

    cpr.setParameters(jesParameters.map(_.toGoogleParamter).toVector.asJava)

    println(s"Pipeline parameters are ${cpr.getParameters}")
    val pipelineId = jesConnection.genomics.pipelines().create(cpr).execute().getPipelineId
    println(s"Pipeline ID is $pipelineId")
    new Pipeline(command, pipelineId, projectId, gcsPath, call, jesParameters, jesConnection.genomics)
  }
}

// Note that id is the JES id not the workflow id
case class Pipeline(command: String,
                    id: String,
                    projectId: String,
                    gcsPath: String,
                    call: Call,
                    jesParameters: Seq[JesParameter],
                    genomicsService: Genomics) {
  def run: Run = Run(this)
}
