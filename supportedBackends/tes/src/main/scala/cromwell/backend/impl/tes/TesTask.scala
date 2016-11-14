package cromwell.backend.impl.tes

import java.nio.file.{FileSystems, Paths}

import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlSingleFile, WdlValue}
import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, OutputEvaluator}
import wdl4s.parser.MemoryUnit

import scala.util.Try

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor) {

  import TesTask._

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val jobPaths = new JobPaths(
    workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key
  )
  private val callEngineFunction = SharedFileSystemExpressionFunctions(
    jobPaths, List(FileSystems.getDefault)
  )
  private val runtimeAttributes: TesRuntimeAttributes = {
    val lookup = jobDescriptor.inputs.apply _
    val evaluateAttrs = jobDescriptor.call.task
      .runtimeAttributes
      .attrs
      .mapValues(_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    TesRuntimeAttributes(runtimeMap, jobDescriptor.workflowDescriptor.workflowOptions)
  }
  private val tesPaths = new TesPaths(jobPaths, runtimeAttributes)

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
      tesPaths.toContainerPath
    )
    // TODO remove this .get and handle error appropriately
    .get

   val inputs = jobDescriptor
    .inputs
    .toSeq
    .flatMap(flattenWdlValueMap)
    .map {
      case (inputName, f: WdlSingleFile) => TaskParameter(
        inputName,
        None,
        tesPaths.storageInput(f.value),
        tesPaths.toContainerPath(f).toString,
        "file",
        false
      )
    }

  val outputs = OutputEvaluator
    .evaluateOutputs(jobDescriptor, callEngineFunction, mapOutputs)
    // TODO remove this .get and handle error appropriately
    .get
    .toSeq
    .map { case (k, v) => (k, v.wdlValue) }
    .flatMap(flattenWdlValueMap)
    .map {
      // TODO handle globs
      case (outputName, WdlSingleFile(path)) => TaskParameter(
        outputName,
        None,
        tesPaths.storagePath(path),
        tesPaths.containerOutput(path),
        "file",
        false
      )
    }

  val workingDirVolume = runtimeAttributes
    .dockerWorkingDir
    .map(path => Volume(
      path,
      // TODO all volumes currently get the same requirements
      Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
      None,
      path
    ))

  val volumes = Seq(
    Volume(
      tesPaths.containerWorkflowRoot,
      Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
      None,
      tesPaths.containerWorkflowRoot
    )
  ) ++ workingDirVolume

  // TODO resolve TES schema around memory format Int -> Double
  val resources = Resources(
    runtimeAttributes.cpu,
    runtimeAttributes.memory.to(MemoryUnit.GB).amount.toInt,
    false,
    Some(volumes),
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage.get,
    // TODO command shouldn't be wrapped in a subshell
    Seq("/bin/bash", "-c", command),
    runtimeAttributes.dockerWorkingDir,
    tesPaths.containerExec("stdout"),
    tesPaths.containerExec("stderr"),
    None
  ))
}

object TesTask {

  private def flattenWdlValueMap(pair: (String, WdlValue)): Seq[(String, WdlValue)] = {
    pair match {
      case (name, file: WdlFile) => Seq((name, file))
      case (name, array: WdlArray) => array.value.zipWithIndex.flatMap {
        case (v: WdlValue, i: Int) => flattenWdlValueMap((name + "-" + i, v))
      }
      case (name, map: WdlMap) => {
        map.value.toSeq flatMap {
          case (name: WdlValue, item: WdlValue) => {
            flattenWdlValueMap((name + name.valueString, item))
          }
        }
      }
      case (name, wdlValue) => Seq((name, wdlValue))
    }
  }

  private final class TesPaths(jobPaths: JobPaths, runtimeAttributes: TesRuntimeAttributes) {

    // Utility for converting a WdlValue so that the path is localized to the
    // container's filesystem.
    def toContainerPath(path: WdlValue): WdlValue = {
      path match {
        case file: WdlFile => {
          val localPath = Paths.get(file.valueString).toAbsolutePath
          WdlFile(containerInput(localPath.toString))
        }
        case array: WdlArray => WdlArray(array.wdlType, array.value map toContainerPath)
        case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toContainerPath)
        case wdlValue => wdlValue
      }
    }

    private def prefixScheme(path: String): String = "file://" + path

    def storageInput(path: String): String = prefixScheme(path)

    // Given an output path, return a path localized to the storage file system
    def storagePath(path: String): String = {
      prefixScheme(jobPaths.callExecutionRoot.resolve(path).toString)
    }

    def containerInput(path: String): String = {
      jobPaths.callDockerRoot.resolve("inputs").resolve(path).toString
    }

    // Given an output path, return a path localized to the container file system
    def containerOutput(path: String): String = containerExec(path)
    // TODO this could be used to create a separate directory for outputs e.g.
    // callDockerRoot.resolve("outputs").resolve(name).toString

    // Given an file name, return a path localized to the container's execution directory
    def containerExec(name: String): String = runtimeAttributes.dockerWorkingDir match {
      case Some(path) => Paths.get(path).resolve(name).toString
      case None => jobPaths.callExecutionDockerRoot.resolve(name).toString
    }

    // The path to the workflow root directory, localized to the container's file system
    val containerWorkflowRoot = jobPaths.dockerWorkflowRoot.toString
  }

  // Utility for converting a WdValue representing an output file path to a WdlValue with
  // a path localized to ____?
  // TODO this is a placeholder for now, until I can fill in the blank
  private def mapOutputs(value: WdlValue) = Try {
    value
  }
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
