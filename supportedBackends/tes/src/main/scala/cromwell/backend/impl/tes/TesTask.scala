package cromwell.backend.impl.tes

import java.nio.file.{FileSystems, Paths}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlSingleFile, WdlValue}
import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, OutputEvaluator}
import cromwell.core.JobOutput
import wdl4s.parser.MemoryUnit
import scala.util.Try

object TesTask {
  private def toDockerPath(path: WdlValue)(implicit jobPaths: JobPaths): WdlValue = {
    path match {
      case file: WdlFile => {
        val localPath = Paths.get(file.valueString).toAbsolutePath
        WdlFile(
          jobPaths.callDockerRoot
          .resolve("inputs")
          .resolve(localPath)
          .toString
        )
      }
      case array: WdlArray => WdlArray(array.wdlType, array.value map toDockerPath)
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toDockerPath)
      case wdlValue => wdlValue
    }
  }

  private def mapOutputs(value: WdlValue) = Try {
    value
  }
}

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor) {

  import TesTask._

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

  // TODO name should be fully qualified
  val name = jobDescriptor.call.fullyQualifiedName
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
    // TODO remove this .get and handle error appropriately
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

  // TODO make paths "file://" URIs

  val outputs = OutputEvaluator
    .evaluateOutputs(jobDescriptor, callEngineFunction, mapOutputs)
    // TODO remove this .get and handle error appropriately
    .get
    .toSeq
    .map {
      // TODO why doesn't WdlValue have unapply
      case (outputName, JobOutput(wdlValue)) => TaskParameter(
        outputName,
        None,
        // TODO probably wrong. should be client-side path (aka storage path)
        wdlValue.valueString,
        // TODO wrong should be worker-side path (aka container path)
        "",
        "file",
        false
      )
    }

  val workingDirVolume = runtimeAttributes
    .dockerWorkingDir
    .map(path => Volume(
      Some(path),
      Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
      None,
      Some(path)
    ))

  val volumes = Seq(
    Volume(
      // Name
      jobPaths.DockerRoot.toString,
      // Size in GB
      Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
      // Source
      None,
      // Mount point
      jobPaths.DockerRoot.toString
    )
  ) ++ workingDirVolume

  // TODO - resolve TES schema around memory format Int -> Double
  val resources = Resources(
    // Minimum CPU cores
    runtimeAttributes.cpu,
    // Minimum RAM in GB
    runtimeAttributes.memory.to(MemoryUnit.GB).amount.toInt,
    // Preemptible?
    false,
    Some(volumes),
    // Zones
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage.get,
    // TODO command shouldn't be wrapped in a subshell
    Seq("/bin/bash", "-c", command),
    runtimeAttributes.dockerWorkingDir,
    jobPaths.callExecutionDockerRoot.resolve("stdout").toString,
    jobPaths.callExecutionDockerRoot.resolve("stderr").toString,
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
                                stdout: String,
                                stderr: String,
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

final case class Volume(name: String,
                        sizeGb: Option[Int],
                        source: Option[String],
                        mountPoint: String)
