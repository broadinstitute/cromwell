package cromwell.engine.backend.jes

import cromwell.binding.{Call, WorkflowDescriptor, CallInputs}
import cromwell.binding.values.WdlFile
import scala.collection.JavaConverters._
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.{Parameter, CreatePipelineRequest}
import JesBackend._

object Pipeline {
  def apply(command: String, workflow: WorkflowDescriptor, call: Call, backendInputs: CallInputs, projectId: String, genomicsService: Genomics): Pipeline = {
    println(s"Command line is $command")
    val runtimeInfo = JesRuntimeInfo(command, call)

    val gcsPath = workflow.callDir(call)

    // We'll need to localize any WdlFiles
    val wdlFileInputs: Map[String, String] = backendInputs collect {case (k, v) if v.isInstanceOf[WdlFile] => k -> v.asInstanceOf[WdlFile].asString}
    
    val cpr = new CreatePipelineRequest
    cpr.setProjectId(projectId)
    cpr.setDocker(runtimeInfo.docker)
    cpr.setResources(runtimeInfo.resources) // FIXME: These are still F'd up if you dig into it
    cpr.setName(workflow.name)

    cpr.setParameters(buildParameters(wdlFileInputs).toVector.asJava)

    println(s"Pipeline parameters are ${cpr.getParameters}")
    val pipelineId = genomicsService.pipelines().create(cpr).execute().getPipelineId
    println(s"Pipeline ID is $pipelineId")
    new Pipeline(command, pipelineId, projectId, gcsPath, call, wdlFileInputs, genomicsService)
  }

  // Combine the stdout/stderr with the job's inputs
  def buildParameters(inputs: Map[String, String]): Iterable[Parameter] = {
     val z =  inputs map {case (name, value) => new Parameter().setName(name).setValue(value).setType("REFERENCE")}
val x = z.toVector
    StdParameters ++ x
  }

  // FIXME: For now we want to always redirect stdout and stderr. This could be problematic if that's what the WDL calls stuff, but oh well
  val StdParameters = Vector(
    new Parameter().setName("stdout").setValue(LocalStdout).setType("REFERENCE"),
    new Parameter().setName("stderr").setValue(LocalStderr).setType("REFERENCE")
  )
}

// FIXME: Members will probably need to change
// Note that id is the JES id not the workflow id
case class Pipeline(command: String,
                    id: String,
                    projectId: String,
                    gcsPath: String,
                    call: Call,
                    inputsToLocalize: Map[String, String],
                    genomicsService: Genomics) {
  def run: Run = Run(this)
}
