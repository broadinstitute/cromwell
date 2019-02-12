package cromwell.backend.impl.tes

import common.collections.EnhancedCollections._
import common.util.StringUtil._
import cromwell.backend.impl.tes.OutputMode.OutputMode
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, Path}
import wdl.draft2.model.FullyQualifiedName
import wdl4s.parser.MemoryUnit
import wom.InstantiatedCommand
import wom.callable.Callable.OutputDefinition
import wom.expression.NoIoFunctionSet
import wom.values._

import scala.language.postfixOps
import scala.util.Try

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor,
                         jobLogger: JobLogger,
                         tesPaths: TesJobPaths,
                         runtimeAttributes: TesRuntimeAttributes,
                         containerWorkDir: Path,
                         commandScriptContents: String,
                         instantiatedCommand: InstantiatedCommand,
                         dockerImageUsed: String,
                         mapCommandLineWomFile: WomFile => WomFile,
                         jobShell: String,
                         outputMode: OutputMode) {

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val workflowName = workflowDescriptor.callable.name
  private val fullyQualifiedTaskName = jobDescriptor.taskCall.fullyQualifiedName
  val name: String = fullyQualifiedTaskName
  val description: String = jobDescriptor.toString

  // TODO validate "project" field of workflowOptions
  val project = {
    workflowDescriptor.workflowOptions.getOrElse("project", "")
  }

  // contains the script to be executed
  private val commandScript = Input(
    name = Option("commandScript"),
    description = Option(fullyQualifiedTaskName + ".commandScript"),
    url = Option(tesPaths.script.pathAsString),
    path = tesPaths.callExecutionDockerRoot.resolve("script").toString,
    `type` = Option("FILE"),
    content = None
  )

  private val commandScriptOut = Output(
    name = Option("commandScript"),
    description = Option(fullyQualifiedTaskName + ".commandScript"),
    url = Option(tesPaths.script.toString),
    path = tesPaths.callExecutionDockerRoot.resolve("script").toString,
    `type` = Option("FILE")
  )

  private def writeFunctionFiles: Map[FullyQualifiedName, Seq[WomFile]] =
    instantiatedCommand.createdFiles map { f => f.file.value.md5SumShort -> List(f.file) } toMap

  private val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor
    .fullyQualifiedInputs
    .safeMapValues {
      _.collectAsSeq { case w: WomFile => w }
    }

  lazy val inputs: Seq[Input] = {
    val result = (callInputFiles ++ writeFunctionFiles).flatMap {
      case (fullyQualifiedName, files) => files.flatMap(_.flattenFiles).zipWithIndex.map {
        case (f, index) =>
          val inputType = f match {
            case _: WomUnlistedDirectory => "DIRECTORY"
            case _: WomSingleFile => "FILE"
            case _: WomGlobFile => "FILE"
          }
          Input(
            name = Option(fullyQualifiedName + "." + index),
            description = Option(workflowName + "." + fullyQualifiedName + "." + index),
            url = Option(f.value),
            path = mapCommandLineWomFile(f).value,
            `type` = Option(inputType),
            content = None
          )
      }
    }.toList ++ Seq(commandScript)
    jobLogger.info(s"Calculated TES inputs (found ${result.size}): " + result.mkString(System.lineSeparator(),System.lineSeparator(),System.lineSeparator()))
    result
  }

  // TODO add TES logs to standard outputs
  private lazy val standardOutputs = Seq("rc", "stdout", "stderr").map {
    f =>
      Output(
        name = Option(f),
        description = Option(fullyQualifiedTaskName + "." + f),
        url = Option(tesPaths.storageOutput(f)),
        path = tesPaths.containerOutput(containerWorkDir, f),
        `type` = Option("FILE")
      )
  }

  // TODO extract output file variable names and match with Files below
  // The problem is that we only care about the files CREATED, so stdout and input redirects are ignored and
  // thus we can't directly match the names returned here to the files returned below. Also we have to consider Arrays
  //
  //  private val outputFileNames = jobDescriptor.call.task.outputs
  //    .filter(o => o.womType.toWdlString == "Array[File]" || o.womType.toWdlString == "File")
  //    .map(_.unqualifiedName)

  // extract output files
  // if output paths are absolute we will ignore them here and assume they are redirects
  private val outputWomFiles: Seq[WomFile] = {
    import cats.syntax.validated._
    // TODO WOM: this should be pushed back into WOM.
    // It's also a mess, evaluateFiles returns an ErrorOr but can still throw. We might want to use an EitherT, although
    // if it fails we just want to fallback to an empty list anyway...
    def evaluateFiles(output: OutputDefinition): List[WomFile] = {
      Try (
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList map { _.file })
      ).getOrElse(List.empty[WomFile].validNel)
       .getOrElse(List.empty)
    }
    
    jobDescriptor.taskCall.callable.outputs
      .flatMap(evaluateFiles)
      .filter(o => !DefaultPathBuilder.get(o.valueString).isAbsolute)
  }

  def handleGlobFile(g: WomGlobFile, index: Int) = {
    val globName = GlobFunctions.globName(g.value)
    val globDirName = "globDir." + index
    val globDirectory = globName + "/"
    val globListName =  "globList." + index
    val globListFile = globName + ".list"
    Seq(
      Output(
        name = Option(globDirName),
        description = Option(fullyQualifiedTaskName + "." + globDirName),
        url = Option(tesPaths.storageOutput(globDirectory)),
        path = tesPaths.containerOutput(containerWorkDir, globDirectory),
        `type` = Option("DIRECTORY")
      ),
      Output(
        name  = Option(globListName),
        description = Option(fullyQualifiedTaskName + "." + globListName),
        url = Option(tesPaths.storageOutput(globListFile)),
        path = tesPaths.containerOutput(containerWorkDir, globListFile),
        `type` = Option("FILE")
      )
    )
  }

  private val womOutputs = outputWomFiles.flatMap(_.flattenFiles)
    .zipWithIndex
    .flatMap {
      case (f: WomSingleFile, index) =>
        val outputFile = f.value
        Seq(
          Output(
            name = Option(fullyQualifiedTaskName + ".output." + index),
            description = Option(fullyQualifiedTaskName + ".output." + index),
            url = Option(tesPaths.storageOutput(outputFile)),
            path = tesPaths.containerOutput(containerWorkDir, outputFile),
            `type` = Option("FILE")
          )
        )
      case (g: WomGlobFile, index) => handleGlobFile(g, index)
      case (d: WomUnlistedDirectory, index) =>
        val directoryPathName = "dirPath." + index
        val directoryPath = d.value.ensureSlashed
        val directoryListName =  "dirList." + index
        val directoryList = d.value.ensureUnslashed + ".list"
        Seq(
          Output(
            name = Option(directoryPathName),
            description = Option(fullyQualifiedTaskName + "." + directoryPathName),
            url = Option(tesPaths.storageOutput(directoryPath)),
            path = tesPaths.containerOutput(containerWorkDir, directoryPath),
            `type` = Option("DIRECTORY")
          ),
          Output(
            name  = Option(directoryListName),
            description = Option(fullyQualifiedTaskName + "." + directoryListName),
            url = Option(tesPaths.storageOutput(directoryList)),
            path = tesPaths.containerOutput(containerWorkDir, directoryList),
            `type` = Option("FILE")
          )
        )
    }

  private val additionalGlobOutput = jobDescriptor.taskCall.callable.additionalGlob.toList.flatMap(handleGlobFile(_, womOutputs.size))
  
  private lazy val cwdOutput = Output(
    name = Option("execution.dir.output"),
    description = Option(fullyQualifiedTaskName + "." + "execution.dir.output"),
    url = Option(tesPaths.callExecutionRoot.pathAsString),
    path = containerWorkDir.pathAsString,
    `type` = Option("DIRECTORY")
  )

  val outputs: Seq[Output] = {
    val result =  outputMode match {
      case OutputMode.GRANULAR => standardOutputs ++ Seq(commandScriptOut) ++ womOutputs ++ additionalGlobOutput
      case OutputMode.ROOT => List(cwdOutput) ++ additionalGlobOutput
    }

    jobLogger.info(s"Calculated TES outputs (found ${result.size}): " + result.mkString(System.lineSeparator(),System.lineSeparator(),System.lineSeparator()))

    result
  }

  private val disk :: ram :: _ = Seq(runtimeAttributes.disk, runtimeAttributes.memory) map {
    case Some(x) =>
      Option(x.to(MemoryUnit.GB).amount)
    case None =>
      None
  }

  val resources = Resources(
    cpu_cores = runtimeAttributes.cpu.map(_.value),
    ram_gb = ram,
    disk_gb = disk,
    preemptible = Option(false),
    zones = None
  )

  val executors = Seq(Executor(
    image = dockerImageUsed,
    command = Seq(jobShell, commandScript.path),
    workdir = runtimeAttributes.dockerWorkingDir,
    stdout = Option(tesPaths.containerOutput(containerWorkDir, "stdout")),
    stderr = Option(tesPaths.containerOutput(containerWorkDir, "stderr")),
    stdin = None,
    env = None
  ))
}

// Field requirements in classes below based off GA4GH schema
final case class Task(id: Option[String],
                      state: Option[String],
                      name: Option[String],
                      description: Option[String],
                      inputs: Option[Seq[Input]],
                      outputs: Option[Seq[Output]],
                      resources: Option[Resources],
                      executors: Seq[Executor],
                      volumes: Option[Seq[String]],
                      tags: Option[Map[String, String]],
                      logs: Option[Seq[TaskLog]])

final case class Executor(image: String,
                          command: Seq[String],
                          workdir: Option[String],
                          stdout: Option[String],
                          stderr: Option[String],
                          stdin: Option[String],
                          env: Option[Map[String, String]])

final case class Input(name: Option[String],
                       description: Option[String],
                       url: Option[String],
                       path: String,
                       `type`: Option[String],
                       content: Option[String])

final case class Output(name: Option[String],
                        description: Option[String],
                        url: Option[String],
                        path: String,
                        `type`: Option[String])

final case class Resources(cpu_cores: Option[Int],
                           ram_gb: Option[Double],
                           disk_gb: Option[Double],
                           preemptible: Option[Boolean],
                           zones: Option[Seq[String]])

final case class OutputFileLog(url: String,
                               path: String,
                               size_bytes: Int)

final case class TaskLog(start_time: Option[String],
                         end_time: Option[String],
                         metadata: Option[Map[String, String]],
                         logs: Option[Seq[ExecutorLog]],
                         outputs: Option[Seq[OutputFileLog]],
                         system_logs: Option[Seq[String]])

final case class ExecutorLog(start_time: Option[String],
                             end_time: Option[String],
                             stdout: Option[String],
                             stderr: Option[String],
                             exit_code: Option[Int])

