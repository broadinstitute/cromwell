package cromwell.backend.impl.tes

import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.core.logging.JobLogger
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.path.Obsolete._
import wdl4s.expression.NoFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.values.{WdlArray, WdlFile, WdlGlobFile, WdlMap, WdlSingleFile, WdlValue}

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor,
                         jobLogger: JobLogger,
                         tesPaths: TesJobPaths,
                         runtimeAttributes: TesRuntimeAttributes) {

  import TesTask._

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val pathBuilders = List(DefaultPathBuilder)
  private val callEngineFunction = SharedFileSystemExpressionFunctions(tesPaths, pathBuilders)
  private val workflowName = workflowDescriptor.workflow.unqualifiedName
  private val fullyQualifiedTaskName = jobDescriptor.call.fullyQualifiedName
  val name = fullyQualifiedTaskName
  val description = jobDescriptor.toString

  runtimeAttributes.dockerWorkingDir match {
          case Some(path: String) => tesPaths.containerWorkingDir = Paths.get(path)
          case _ => tesPaths.containerWorkingDir = tesPaths.callExecutionDockerRoot
  }

  // TODO validate "project" field of workflowOptions
  val project = {
    val workflowName = jobDescriptor.workflowDescriptor.rootWorkflow.unqualifiedName
    workflowDescriptor.workflowOptions.getOrElse("project", workflowName)
  }

  private val commandScript = TaskParameter(
    "commandScript",
    Some(fullyQualifiedTaskName + ".commandScript"),
    tesPaths.storageInput(tesPaths.script.toString),
    tesPaths.callExecutionDockerRoot.resolve("script").toString,
    "File",
    Some(false)
  )

  val inputs: Seq[TaskParameter] = jobDescriptor
    .fullyQualifiedInputs
    .toSeq
    .flatMap(flattenWdlValueMap)
    .filter{
        case (_: String, _: WdlSingleFile) => true
        case _ => false
    }
    .map {
      case (inputName, f) => TaskParameter(
        inputName,
        Some(workflowName + "." + inputName),
        tesPaths.storageInput(f.valueString),
        tesPaths.toContainerPath(f).valueString,
        "File",
        Some(false)
      )
    } ++ Seq(commandScript)

  // TODO add TES logs to standard outputs
  private val standardOutputs = Seq("rc", "stdout", "stderr").map {
    f =>
      TaskParameter(
        f,
        Some(fullyQualifiedTaskName + "." + f),
        tesPaths.storageOutput(f),
        tesPaths.containerOutput(f),
        "File",
        Some(false)
      )
  }

  // extract output file variable names
  private val outputFileNames = jobDescriptor.call.task.outputs
      .filter {
        case o => o.wdlType.toWdlString == "Array[File]" || o.wdlType.toWdlString == "File"
      }
    .map{case o => o.unqualifiedName}

  // extract output files
  private val outputWdlFiles = jobDescriptor.call.task
    .findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions)

  private val wdlOutputs = outputFileNames.zip(outputWdlFiles)
    .flatMap {
      case (outputName, f: WdlSingleFile) => {
        val outputFile = f.valueString
        Seq(
          TaskParameter(
            outputName,
            Some(fullyQualifiedTaskName + "." + outputName),
            tesPaths.storageOutput(outputFile),
            tesPaths.containerOutput(outputFile),
            "File",
            Some(false)
          )
        )
      }
      case (outputName, g: WdlGlobFile) => {
        val globName = callEngineFunction.globName(g.value)
        val globDirName = outputName + ".globDir"
        val globDirectory = globName + "/"
        val globListName = outputName + ".globList"
        val globListFile = globName + ".list"
        Seq(
          TaskParameter(
            globDirName,
            Some(fullyQualifiedTaskName + "." + globDirName),
            tesPaths.storageOutput(globDirectory),
            tesPaths.containerOutput(globDirectory),
            "Directory",
            Some(false)
          ),
          TaskParameter(
            globListName,
            Some(fullyQualifiedTaskName + "." + globListName),
            tesPaths.storageOutput(globListFile),
            tesPaths.containerOutput(globListFile),
            "File",
            Some(false)
          )
        )
      }
    }

  val outputs = wdlOutputs ++ standardOutputs

  private val workingDirVolume = runtimeAttributes
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
      tesPaths.dockerWorkflowRoot.toString,
      Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
      None,
      tesPaths.dockerWorkflowRoot.toString
    )
  ) ++ workingDirVolume

  val resources = Resources(
    runtimeAttributes.cpu,
    runtimeAttributes.memory.to(MemoryUnit.GB).amount.toInt,
    Some(false),
    Some(volumes),
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage,
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
      case (_, map: WdlMap) => {
        map.value.toSeq flatMap {
          case (name: WdlValue, item: WdlValue) => {
            flattenWdlValueMap((name + name.valueString, item))
          }
        }
      }
      case (name, wdlValue) => Seq((name, wdlValue))
    }
  }
}


final case class TesTaskMessage(name: Option[String],
                                description: Option[String],
                                projectId: Option[String],
                                inputs: Option[Seq[TaskParameter]],
                                outputs: Option[Seq[TaskParameter]],
                                resources: Resources,
                                docker: Seq[DockerExecutor])

final case class DockerExecutor(imageName: String,
                                cmd: Seq[String],
                                workdir: Option[String],
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
