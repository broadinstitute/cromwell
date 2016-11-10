package cromwell.backend.impl.tes

import java.nio.file.{FileSystems, Paths}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlSingleFile, WdlValue}
import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import wdl4s.parser.MemoryUnit
import scala.util.Try

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
      case file: WdlFile => {
        val localPath = Paths.get(file.valueString)

        localPath.toAbsolutePath match {
          case p if p.startsWith(jobPaths.DockerRoot) => WdlFile(p.toString)
          case p =>
            val fileName = p.getFileName
            val callPath = jobPaths.callRoot.resolve(fileName)
            val subPath = callPath.subpath(jobPaths.executionRoot.getNameCount, callPath.getNameCount)
            WdlFile(jobPaths.DockerRoot.resolve(subPath).toString)
        }
      }
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

  val command = jobDescriptor
    .call
    .instantiateCommandLine(
      jobDescriptor.inputs,
      callEngineFunction,
      toDockerPath
    )
    .get

  val inputs = jobDescriptor
    .inputs
    .toSeq
    .map {
      case (inputName, f: WdlSingleFile) => TaskParameter(
        inputName,
        None,
        f.value,
        toDockerPath(f).toString,
        "file",
        false
      )
    }

  val outputs = Seq()

  val volumes = Some(Seq(
    Volume(
      // Name
      Some("cromwell_inputs"),
      // Size in GB
      Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
      // Source
      None,
      // Mount point
      Some("/tmp")
    )
  ))

  // TODO - resolve TES schema around memory format Int -> Double
  val resources = Resources(
    // Minimum CPU cores
    runtimeAttributes.cpu,
    // Minimum RAM in GB
    runtimeAttributes.memory.to(MemoryUnit.GB).amount.toInt,
    // Preemptible?
    false,
    volumes,
    // Zones
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage.get,
    // TODO command shouldn't be wrapped in a subshell
    Seq("/bin/bash", "-c", command),
    runtimeAttributes.dockerWorkingDir,
    Some("/tmp/stdout"),
    Some("/tmp/stderr"),
    None
  ))
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

final case class Resources(minimumCpuCores: Int,
                           minimumRamGb: Int,
                           preemptible: Boolean,
                           volumes: Option[Seq[Volume]],
                           zones: Option[Seq[String]])

final case class Volume(name: Option[String],
                        sizeGb: Option[Int],
                        source: Option[String],
                        mountPoint: Option[String])
