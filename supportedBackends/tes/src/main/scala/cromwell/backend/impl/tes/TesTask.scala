package cromwell.backend.impl.tes

import java.nio.file.Paths
import better.files.File
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlSingleFile, WdlValue}
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.backend.wdl.OutputEvaluator
import cromwell.core.logging.JobLogger
import cromwell.core.path.DefaultPathBuilder
import lenthall.util.TryUtil
import wdl4s.parser.MemoryUnit
import scala.util.Try

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor,
                         jobLogger: JobLogger) {

  import TesTask._

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val pathBuilders = List(DefaultPathBuilder)
  private val jobPaths = new JobPathsWithDocker(jobDescriptor.key, workflowDescriptor, configurationDescriptor.backendConfig)
  private val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, pathBuilders)

  private val runtimeAttributes = {
    val lookup = jobDescriptor.fullyQualifiedInputs.apply _
    val evaluateAttrs = jobDescriptor.key.call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    TesRuntimeAttributes(runtimeMap, jobDescriptor.workflowDescriptor.workflowOptions)
  }

  private val tesPaths = new TesPaths(jobPaths, runtimeAttributes)

  private val workflowName = workflowDescriptor.workflow.unqualifiedName
  val name = jobDescriptor.call.fullyQualifiedName
  val description = jobDescriptor.toString
  val taskId = jobDescriptor.toString

  // TODO validate "project" field of workflowOptions
  val project = {
    val workflowName = jobDescriptor.workflowDescriptor.rootWorkflow.unqualifiedName
    workflowDescriptor.workflowOptions.getOrElse("project", workflowName)
  }

  val commandString = jobDescriptor
    .key
    .call
    .task
    .instantiateCommand(
      jobDescriptor.inputDeclarations,
      callEngineFunction,
      tesPaths.toContainerPath
    )
    // TODO remove this .get and handle error appropriately
    .get

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(instantiatedCommand: String) = {
    val rcPath = jobPaths.toDockerPath(jobPaths.returnCode)
    val rcTmpPath = s"$rcPath.tmp"

    val scriptBody = s"""

#!/bin/sh
(
 $instantiatedCommand
)
echo $$? > $rcTmpPath
mv $rcTmpPath $rcPath

""".trim + "\n"

    File(jobPaths.script).write(scriptBody)
  }

  File(jobPaths.callExecutionRoot).createDirectories()
  jobLogger.info(s"`\n$commandString`")
  writeScript(commandString)

  private val commandScript = TaskParameter(
    name + ".commandScript",
    None,
    tesPaths.storageInput(jobPaths.script.toString),
    jobPaths.callExecutionDockerRoot.resolve("script").toString,
    "File",
    Some(false)
  )

  val inputs = jobDescriptor
    .fullyQualifiedInputs
    .toSeq
    .flatMap(flattenWdlValueMap)
    .map {
      case (inputName, f: WdlSingleFile) => TaskParameter(
        workflowName + "." + inputName,
        None,
        tesPaths.storageInput(f.valueString),
        tesPaths.toContainerPath(f).valueString,
        "File",
        Some(false)
      )
    } ++ Seq(commandScript)

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
        tesPaths.storageOutput(path),
        tesPaths.containerOutput(path),
        "File",
        Some(false)
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
    Some(false),
    Some(volumes),
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage.get,
    // TODO command shouldn't be wrapped in a subshell
    Seq("/bin/bash", commandScript.path),
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

  private final class TesPaths(jobPaths: JobPathsWithDocker, runtimeAttributes: TesRuntimeAttributes) {

    // Utility for converting a WdlValue so that the path is localized to the
    // container's filesystem.
    def toContainerPath(path: WdlValue): WdlValue = {
      path match {
        case file: WdlFile => {
          val localPath = Paths.get(file.valueString).toAbsolutePath
          val containerPath = containerInput(localPath.toString)
          WdlFile(containerPath)
        }
        case array: WdlArray => WdlArray(array.wdlType, array.value map toContainerPath)
        case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toContainerPath)
        case wdlValue => wdlValue
      }
    }

    private def prefixScheme(path: String): String = "file://" + path

    def storageInput(path: String): String = prefixScheme(path)

    // Given an output path, return a path localized to the storage file system
    def storageOutput(path: String): String = {
      prefixScheme(jobPaths.callExecutionRoot.resolve(path).toString)
    }

    def containerInput(path: String): String = {
      jobPaths.callDockerRoot.resolve("inputs").toString + path
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


final case class TesTaskMessage(name: Option[String],
                                description: Option[String],
                                projectId: Option[String],
                                taskID: String,
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
                               create: Option[Boolean])

final case class Resources(minimumCpuCores: Int,
                           minimumRamGb: Int,
                           preemptible: Option[Boolean],
                           volumes: Option[Seq[Volume]],
                           zones: Option[Seq[String]])

final case class Volume(name: String,
                        sizeGb: Option[Int],
                        source: Option[String],
                        mountPoint: String)
