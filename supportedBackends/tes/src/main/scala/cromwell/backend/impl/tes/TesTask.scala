package cromwell.backend.impl.tes

import java.nio.file.{FileSystems, Paths}
import scala.util.{Failure, Try}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlSingleFile, WdlValue}
import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}


final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor) {

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)
  private val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, List(FileSystems.getDefault))

  private val runtimeAttributes = {
    val lookup = jobDescriptor.inputs.apply _
    val evaluateAttrs = jobDescriptor.call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    TesRuntimeAttributes(runtimeMap, jobDescriptor.workflowDescriptor.workflowOptions)
  }

  private def toDockerPath(path: WdlValue): WdlValue = {
    path match {
      case file: WdlFile => WdlFile(jobPaths.toDockerPath(Paths.get(path.valueString)).toString)
      case array: WdlArray => WdlArray(array.wdlType, array.value map toDockerPath)
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toDockerPath)
      case wdlValue => wdlValue
    }
  }


  val name = jobDescriptor.call.task.name
  val description = jobDescriptor.toString
  val taskId = jobDescriptor.toString

  // TODO validate "project" field of workflowOptions
  val project = {
    val workflowName = jobDescriptor.call.rootWorkflow.unqualifiedName
    workflowDescriptor.workflowOptions.getOrElse("project", workflowName)
  }

  val commands = jobDescriptor
    .call
    .instantiateCommandLine(
      jobDescriptor.inputs
      callEngineFunction,
      toDockerFile
    )

  val inputs = jobDescriptor
    .inputs
    .toSeq
    .map {
      case (inputName, f: WdlSingleFile) => TaskParameter(
        inputName,       // Name
        None,            // Description
        f.value,         // Source path
        toDockerPath(f), // Destination path   NOTE: this must match a Volume.mountPoint in the Resources
        "file",          // Type
        false            // Create?
      )
    })

  val outputs = Seq()

  val volumes = Some(Seq(
    Volume(
      // Name
      Some("cromwell_inputs"),
      // Size in GB
      Some(runtimeAttributes.disk.amount.toInt),
      // Source
      None,
      // Mount point
      Some("/tmp")
    )
  ))

  // TODO - resolve TES schema around memory format Int -> Double
  val resources = Resources(
    // Minimum CPU cores
    runtimeAttributes.cpu
    // Minimum RAM in GB
    runtimeAttribute.memory.amount.toInt,
    // Preemptible?
    None, 
    volumes,
    // Zones
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage.get,
    // TODO command shouldn't be wrapped in a subshell
    Seq("/bin/bash", "-c", cmd),
    runtimeAttributes.dockerWorkingDir,
    Some("/tmp/stdout"),
    Some("/tmp/stderr"),
    None
  )
}


final case class TesTaskMessage(name: String,
                                description: String,
                                projectId: String,
                                taskId: String,
                                inputs: Option[Seq[TaskParameter]],
                                outputs: Option[Seq[TaskParameter]],
                                resources: Resources,
                                docker: Seq[DockerExecutor])

final case class DockerExecutor(imageName: String,
                                cmd: Seq[String],
                                workDir: Option[String],
                                stdout: Option[String],
                                stderr: Option[String],
                                stdin: Option[String])

final case class TaskParameter(name: String,
                               description: Option[String],
                               location: String,
                               path: String,
                               `class`: String,
                               create: Boolean)

final case class Resources(minimumCpuCores: Option[Int],
                           minimumRamGb: Option[Int],
                           preemptible: Option[Boolean],
                           volumes: Option[Seq[Volume]],
                           zones: Option[Seq[String]])

final case class Volume(name: Option[String],
                        sizeGb: Option[Int],
                        source: Option[String],
                        mountPoint: Option[String])
