package cromwell.backend.impl.tes

import cromwell.backend.standard.StandardExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, Path}
import wdl4s.FullyQualifiedName
import wdl4s.expression.NoFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.values.{WdlFile, WdlGlobFile, WdlSingleFile, WdlValue}

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor,
                         jobLogger: JobLogger,
                         tesPaths: TesJobPaths,
                         runtimeAttributes: TesRuntimeAttributes,
                         containerWorkDir: Path,
                         backendEngineFunctions: StandardExpressionFunctions) {

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val workflowName = workflowDescriptor.workflow.unqualifiedName
  private val fullyQualifiedTaskName = jobDescriptor.call.fullyQualifiedName
  val name: String = fullyQualifiedTaskName
  val description: String = jobDescriptor.toString

  // TODO validate "project" field of workflowOptions
  val project = {
    workflowDescriptor.workflowOptions.getOrElse("project", workflowName)
  }

  // contains the script to be executed
  private val commandScript = TaskParameter(
    Option("commandScript"),
    Option(fullyQualifiedTaskName + ".commandScript"),
    tesPaths.storageInput(tesPaths.script.toString),
    tesPaths.callExecutionDockerRoot.resolve("script").toString,
    "File",
    Option(false)
  )

  private def writeFunctionFiles(commandLineValueMapper: WdlValue => WdlValue): Map[FullyQualifiedName, Seq[WdlFile]] = {
    val commandLineMappedInputs = jobDescriptor.inputDeclarations map {
      case (declaration, value) => declaration.fullyQualifiedName -> commandLineValueMapper(value)
    }

    jobDescriptor
      .call
      .task
      .evaluateFilesFromCommand(commandLineMappedInputs, backendEngineFunctions)
      .map {
        case (expression, file) => expression.toWdlString -> Seq(file)
      }
  }

  private val callInputFiles: Map[FullyQualifiedName, Seq[WdlFile]] = jobDescriptor
    .fullyQualifiedInputs
    .mapValues {
      _.collectAsSeq { case w: WdlFile => w }
    }

  def inputs(commandLineValueMapper: WdlValue => WdlValue): Seq[TaskParameter] = (callInputFiles ++ writeFunctionFiles(commandLineValueMapper))
    .flatMap {
      case (fullyQualifiedName, files) => files.zipWithIndex.map {
        case (f, index) => TaskParameter(
          Option(fullyQualifiedName + "." + index),
          Option(workflowName + "." + fullyQualifiedName + "." + index),
          tesPaths.storageInput(f.value),
          tesPaths.containerInput(f.value),
          "File",
          Option(false)
        )
      }
    }.toList ++ Seq(commandScript)

  // TODO add TES logs to standard outputs
  private val standardOutputs = Seq("rc", "stdout", "stderr").map {
    f =>
      TaskParameter(
        Option(f),
        Option(fullyQualifiedTaskName + "." + f),
        tesPaths.storageOutput(f),
        tesPaths.containerOutput(containerWorkDir, f),
        "File",
        Option(false)
      )
  }

  // TODO extract output file variable names and match with Files below
  // The problem is that we only care about the files CREATED, so stdout and input redirects are ignored and
  // thus we can't directly match the names returned here to the files returned below. Also we have to consider Arrays
  //
  //  private val outputFileNames = jobDescriptor.call.task.outputs
  //    .filter(o => o.wdlType.toWdlString == "Array[File]" || o.wdlType.toWdlString == "File")
  //    .map(_.unqualifiedName)

  // extract output files
  // if output paths are absolute we will ignore them here and assume they are redirects
  private val outputWdlFiles: Seq[WdlFile] = jobDescriptor.call.task
    .findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions)
    .filter(o => !DefaultPathBuilder.get(o.valueString).isAbsolute)

  private val wdlOutputs = outputWdlFiles
    .zipWithIndex
    .flatMap {
      case (f: WdlSingleFile, index) =>
        val outputFile = f.value
        Seq(
          TaskParameter(
            Option(fullyQualifiedTaskName + ".output." + index),
            Option(fullyQualifiedTaskName + ".output." + index),
            tesPaths.storageOutput(outputFile),
            tesPaths.containerOutput(containerWorkDir, outputFile),
            "File",
            Option(false)
          )
        )
      case (g: WdlGlobFile, index) =>
        val globName = backendEngineFunctions.globName(g.value)
        val globDirName = "globDir." + index
        val globDirectory = globName + "/"
        val globListName =  "globList." + index
        val globListFile = globName + ".list"
        Seq(
          TaskParameter(
            Option(globDirName),
            Option(fullyQualifiedTaskName + "." + globDirName),
            tesPaths.storageOutput(globDirectory),
            tesPaths.containerOutput(containerWorkDir, globDirectory),
            "Directory",
            Option(false)
          ),
          TaskParameter(
            Option(globListName),
            Option(fullyQualifiedTaskName + "." + globListName),
            tesPaths.storageOutput(globListFile),
            tesPaths.containerOutput(containerWorkDir, globListFile),
            "File",
            Option(false)
          )
        )
    }

  val outputs: Seq[TaskParameter] = wdlOutputs ++ standardOutputs

  // TODO all volumes currently get the same disk requirements
  private val workingDirVolume = runtimeAttributes
    .dockerWorkingDir
    .map(path => Volume(
      Option(path),
      runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt,
      None,
      path,
      Option(false)
    ))

  val volumes = Seq(
    Volume(
      Option(tesPaths.callInputsDockerRoot.toString),
      runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt,
      None,
      tesPaths.callInputsDockerRoot.toString,
      // inputs in read-only volume
      Option(true)
    ),
    Volume(
      Option(tesPaths.callExecutionDockerRoot.toString),
      runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt,
      None,
      tesPaths.callExecutionDockerRoot.toString,
      Option(false)
    )
  ) ++ workingDirVolume

  val resources = Resources(
    runtimeAttributes.cpu,
    runtimeAttributes.memory.to(MemoryUnit.GB).amount.toInt,
    Option(false),
    volumes,
    None
  )

  val dockerExecutor = Seq(DockerExecutor(
    runtimeAttributes.dockerImage,
    Seq("/bin/bash", commandScript.path),
    runtimeAttributes.dockerWorkingDir,
    Option(tesPaths.containerOutput(containerWorkDir, "stdout")),
    Option(tesPaths.containerOutput(containerWorkDir, "stderr")),
    None,
    None
  ))
}

// Field requirements in classes below based off GA4GH schema

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
                                stdout: Option[String],
                                stderr: Option[String],
                                stdin: Option[String],
                                ports: Option[Seq[Ports]])

final case class TaskParameter(name: Option[String],
                               description: Option[String],
                               location: String,
                               path: String,
                               `class`: String,
                               create: Option[Boolean])

final case class Resources(minimumCpuCores: Int,
                           minimumRamGb: Int,
                           preemptible: Option[Boolean],
                           volumes: Seq[Volume],
                           zones: Option[Seq[String]])

final case class Volume(name: Option[String],
                        sizeGb: Int,
                        source: Option[String],
                        mountPoint: String,
                        readonly: Option[Boolean])

final case class JobLogs(cmd: Option[Seq[String]],
                         startTime: Option[String],
                         endTime: Option[String],
                         stdout: Option[String],
                         stderr: Option[String],
                         exitCode: Option[Int],
                         hostIP: Option[String],
                         ports: Option[Seq[Ports]])

final case class Ports(host: Option[String],
                       container: String)