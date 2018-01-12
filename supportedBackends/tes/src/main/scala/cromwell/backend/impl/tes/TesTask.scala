package cromwell.backend.impl.tes

import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.core.NoIoFunctionSet
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, Path}
import wdl.FullyQualifiedName
import wdl4s.parser.MemoryUnit
import wom.InstantiatedCommand
import wom.callable.Callable.OutputDefinition
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
                         dockerImageUsed: String) {

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
    url = None,
    path = tesPaths.callExecutionDockerRoot.resolve("script").toString,
    `type` = Option("FILE"),
    content = Option(commandScriptContents)
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
    .mapValues {
      _.collectAsSeq { case w: WomFile => w }
    }

  def inputs(commandLineValueMapper: WomValue => WomValue): Seq[Input] =
    (callInputFiles ++ writeFunctionFiles).flatMap {
      case (fullyQualifiedName, files) => files.zipWithIndex.map {
        case (f, index) => Input(
          name = Option(fullyQualifiedName + "." + index),
          description = Option(workflowName + "." + fullyQualifiedName + "." + index),
          url = Option(f.value),
          path = tesPaths.containerInput(f.value),
          `type` = Option("FILE"),
          content = None
        )
      }
    }.toList ++ Seq(commandScript)

  // TODO add TES logs to standard outputs
  private val standardOutputs = Seq("rc", "stdout", "stderr").map {
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
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList)
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

  private val womOutputs = outputWomFiles
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
      case (unsupported: WomFile, _) =>
        // TODO: WOM: WOMFILE: Add support for directories.
        throw new NotImplementedError(s"$unsupported is not supported yet.")
    }
  
  private val additionalGlobOutput = jobDescriptor.taskCall.callable.additionalGlob.toList.flatMap(handleGlobFile(_, womOutputs.size))

  val outputs: Seq[Output] = womOutputs ++ standardOutputs ++ Seq(commandScriptOut) ++ additionalGlobOutput

  private val disk :: ram :: _ = Seq(runtimeAttributes.disk, runtimeAttributes.memory) map {
    case Some(x) =>
      Option(x.to(MemoryUnit.GB).amount)
    case None =>
      None
  }

  val resources = Resources(
    cpu_cores = runtimeAttributes.cpu,
    ram_gb = ram,
    disk_gb = disk,
    preemptible = Option(false),
    zones = None
  )

  val executors = Seq(Executor(
    image = dockerImageUsed,
    command = Seq("/bin/bash", commandScript.path),
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

