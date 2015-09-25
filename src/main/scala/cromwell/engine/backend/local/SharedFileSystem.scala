package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Files, Path, Paths}
import java.security.MessageDigest

import com.typesafe.config.ConfigFactory
import cromwell.binding._
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlMapType}
import cromwell.binding.values.{WdlValue, _}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.backend.{LocalFileSystemBackendCall, StdoutStderr}
import cromwell.engine.db.DataAccess
import org.apache.commons.io.FileUtils

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystem {
  val SharedFileSystemConf = ConfigFactory.load.getConfig("backend").getConfig("shared-filesystem")
  val CromwellExecutionRoot = SharedFileSystemConf.getString("root")
  val LocalizationStrategies = SharedFileSystemConf.getStringList("localization").asScala
  val Localizers = LocalizationStrategies map {
    case "hard-link" => localizePathViaHardLink _
    case "soft-link" => localizePathViaSymbolicLink _
    case "copy" => localizePathViaCopy _
    case unsupported => throw new UnsupportedOperationException(s"Localization strategy $unsupported is not recognized")
  }

  private def localizePathViaCopy(call: Option[Call], originalPath: Path, executionPath: Path): Try[Unit] = {
    if (Files.isDirectory(originalPath)) {
      Try(FileUtils.copyDirectory(originalPath.toFile, executionPath.toFile))
    } else {
      Try(FileUtils.copyFile(originalPath.toFile, executionPath.toFile))
    }
  }

  private def localizePathViaHardLink(call: Option[Call], originalPath: Path, executionPath: Path): Try[Unit] =
    Try(Files.createLink(executionPath, originalPath))

  /**
   * TODO: The 'call' parameter here represents the call statement in WDL that references this path.
   * We're supposed to not use symbolic links if the call uses Docker.  However, this is currently a
   * bit incorrect because multiple calls can reference the same path if that path is in a declaration.
   *
   * The symbolic link will only fail in the Docker case if a Call uses the file directly and not
   * indirectly through one of its input expressions
   */
  private def localizePathViaSymbolicLink(call: Option[Call], originalPath: Path, executionPath: Path): Try[Unit] = {
    call.flatMap(_.docker) match {
      case Some(_) => Failure(new UnsupportedOperationException("Cannot localize with symbolic links with Docker"))
      case _ => Try(Files.createSymbolicLink(executionPath, originalPath.toAbsolutePath))
    }
  }
}

trait SharedFileSystem {

  import SharedFileSystem._

  def postProcess(backendCall: LocalFileSystemBackendCall): Try[CallOutputs] = {
    // Evaluate output expressions, performing conversions from String -> File where required.
    val outputMappings = backendCall.call.task.outputs map { taskOutput =>
      val tryConvertedValue =
        for {
          expressionValue <- taskOutput.expression.evaluate(backendCall.lookupFunction, backendCall.engineFunctions)
          convertedValue <- outputAutoConversion(backendCall, taskOutput, expressionValue)
          pathAdjustedValue <- Success(absolutizeOutputWdlFile(convertedValue, backendCall.callRootPath))
        } yield pathAdjustedValue
      taskOutput.name -> tryConvertedValue
    }

    val taskOutputFailures = outputMappings filter { _._2.isFailure }

    if (taskOutputFailures.isEmpty) {
      val unwrappedMap = outputMappings collect { case (name, Success(wdlValue)) => name -> wdlValue }
      Success(unwrappedMap.toMap)
    } else {
      val message = taskOutputFailures collect { case (name, Failure(e)) => s"$name: $e" }
      Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: ${message.mkString("\n")}"))
    }
  }

  def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): StdoutStderr = {
    val dir = LocalBackend.hostCallPath(descriptor.namespace.workflow.name, descriptor.id, callName, index)
    StdoutStderr(
      stdout = WdlFile(dir.resolve("stdout").toAbsolutePath.toString),
      stderr = WdlFile(dir.resolve("stderr").toAbsolutePath.toString)
    )
  }

  /**
   * Creates host execution directory, inputs path, and outputs path.  Stages any input files into the workflow-inputs
   * directory and localizes their paths relative to the container.
   */
  def initializeForWorkflow(descriptor: WorkflowDescriptor): Try[HostInputs] = {
    val hostExecutionDirectory = LocalBackend.hostExecutionPath(descriptor).toFile
    hostExecutionDirectory.mkdirs()
    val hostExecutionAbsolutePath = hostExecutionDirectory.getAbsolutePath
    Array("workflow-inputs", "workflow-outputs") foreach { Paths.get(hostExecutionAbsolutePath, _).toFile.mkdir() }
    stageWorkflowInputs(descriptor)
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = {
    def containerPath(path: String): WdlFile = {
      // Host path would look like cromwell-executions/three-step/f00ba4/call-ps/stdout.txt
      // Container path should look like /root/f00ba4/call-ps/stdout.txt
      val fullPath = Paths.get(path).toFile.getAbsolutePath
      // Strip out everything before cromwell-executions.
      val pathUnderCromwellExecutions = fullPath.substring(fullPath.indexOf(CromwellExecutionRoot) + CromwellExecutionRoot.length)
      // Strip out the workflow name (the first component under cromwell-executions).
      val pathWithWorkflowName = Paths.get(pathUnderCromwellExecutions)
      WdlFile(WdlFile.appendPathsWithSlashSeparators("/root", pathWithWorkflowName.subpath(1, pathWithWorkflowName.getNameCount).toString))
    }

    // If this call is using Docker, adjust input paths, otherwise return unaltered input paths.
    def adjustPath(nameAndValue: (String, WdlValue)): (String, WdlValue) = {
      val (name, value) = nameAndValue
      val adjusted = value match {
        case WdlFile(path) => containerPath(path)
        case WdlArray(t, values) => new WdlArray(t, values map { adjustPath(name, _)._2 })
        case WdlMap(t, values) => new WdlMap(t, values mapValues { adjustPath(name, _)._2 })
        case x => x
      }
      name -> adjusted
    }

    // If this call is using Docker, adjust input paths, otherwise return unaltered input paths.
    if (call.docker.isDefined) inputs map adjustPath else inputs
  }

  /**
   * Given the specified workflow descriptor and inputs, stage any WdlFiles into the workflow-inputs subdirectory
   * of the workflow execution directory.  Return a Map of the input values with any input WdlFiles adjusted to
   * reflect host paths.
   *
   * Original input path: /could/be/anywhere/input.bam
   * Host inputs path: $PWD/cromwell-executions/some-workflow-name/0f00-ba4/workflow-inputs/input.bam
   */
  private def stageWorkflowInputs(descriptor: WorkflowDescriptor): Try[HostInputs] = {
    val hostInputsPath = Paths.get(LocalBackend.hostExecutionPath(descriptor).toFile.getAbsolutePath, "workflow-inputs")
    val attemptedStagedInputs = descriptor.actualInputs map { case(name, value) =>
      val call = descriptor.namespace.resolve(name.split("\\.").dropRight(1).mkString(".")) match {
        case Some(c: Call) => Some(c)
        case _ => None
      }
      name -> stageWdlValue(call, value, hostInputsPath)
    }
    attemptedStagedInputs.values collect { case f: Failure[_] => f } match {
      case Nil =>
        Success(attemptedStagedInputs map { case (k, v) => k -> v.get })
      case failedInputs =>
        val errors = failedInputs map { _.exception.getMessage } mkString "\n"
        Failure(new UnsupportedOperationException(s"Failures during localization:\n\n$errors"))
    }
  }

  private def stageInputArray(call: Option[Call], array: WdlArray, hostInputsPath: Path): Try[WdlArray] = {
    val attemptedStagedArray = array.value.map(stageWdlValue(call, _, hostInputsPath))
    attemptedStagedArray collect { case f: Failure[_] => f } match {
      case s: Seq[Failure[_]] if s.nonEmpty =>
        val errors = s map { _.exception.getMessage } mkString "\n"
        Failure(new UnsupportedOperationException(s"Failures during localization of array $array:\n\n$errors"))
      case _ => Success(WdlArray(array.wdlType, attemptedStagedArray.map(_.get)))
    }
  }

  private def stageWdlValue(call: Option[Call], value: WdlValue, hostInputsPath: Path): Try[WdlValue] = value match {
    case w: WdlFile => stageWdlFile(call, w, hostInputsPath)
    case a: WdlArray => stageInputArray(call, a, hostInputsPath)
    case x => Success(x)
  }

  private def stageWdlFile(call: Option[Call], wdlFile: WdlFile, hostInputsPath: Path): Try[WdlFile] = {
    val originalPath = Paths.get(wdlFile.value)
    val directoryIdentifier = MessageDigest.getInstance("MD5").digest(originalPath.toAbsolutePath.getParent.toString.getBytes) map {byte => f"$byte%02x"} mkString
    val executionPath = hostInputsPath.resolve(s"${directoryIdentifier.substring(0,8)}-${originalPath.getFileName.toString}")

    val attemptedLocalization = Stream(Localizers: _*) map { _(call, originalPath, executionPath) } find { _.isSuccess }
    attemptedLocalization match {
      case Some(_) => Success(WdlFile(executionPath.toString))
      case None => Failure(throw new UnsupportedOperationException(s"Could not localize $wdlFile -> $executionPath"))
    }
  }

  private def assertTaskOutputPathExists(path: String, taskOutput: TaskOutput, callFqn: String): Try[WdlFile] =
    if (Files.exists(Paths.get(path))) Success(WdlFile(path))
    else Failure(new RuntimeException(
      s"""ERROR: Could not process output '${taskOutput.wdlType.toWdlString} ${taskOutput.name}' of $callFqn:
         |
         |Invalid path: $path
       """.stripMargin
    ))

  private def hostAbsoluteFilePath(backendCall: LocalFileSystemBackendCall, pathString: String): String =
    if (new File(pathString).isAbsolute) pathString else Paths.get(backendCall.callRootPath.toAbsolutePath.toString, pathString).toString

  /**
   * Handle possible auto-conversion from an output expression `WdlString` to a `WdlFile` task output.
   * The following should work:
   *
   * <pre>
   * File bam = "foo.bam"
   * </pre>
   *
   * This also handles coercions to target types.  For example, a task output may look like this:
   *
   * <pre>
   * Map[Int, String] my_map = read_map(stdout())
   * </pre>
   *
   * read_map() will return a Map[String, String] but if the target type is Map[Int, String], this
   * function will attempt to do the coercion.
   *
   * read_lines() will return an Array[String] but if the target type is Array[Float], then this
   * function will do that conversion.
   */
  private def outputAutoConversion(backendCall: LocalFileSystemBackendCall, taskOutput: TaskOutput, rawOutputValue: WdlValue): Try[WdlValue] = {
    rawOutputValue match {
      case v: WdlString if taskOutput.wdlType == WdlFileType => assertTaskOutputPathExists(hostAbsoluteFilePath(backendCall, v.value), taskOutput, backendCall.call.fullyQualifiedName)
      case m: WdlMap if taskOutput.wdlType.isInstanceOf[WdlMapType] => taskOutput.wdlType.coerceRawValue(m)
      case a: WdlArray if taskOutput.wdlType.isInstanceOf[WdlArrayType] => taskOutput.wdlType.coerceRawValue(a)
      case v if v.wdlType == taskOutput.wdlType => Success(v)
      case _ => Failure(new RuntimeException(
          s"""Error processing '${backendCall.call.fullyQualifiedName}.${taskOutput.name}':
             |
             |Value $rawOutputValue cannot be converted to ${taskOutput.wdlType.toWdlString}
           """.stripMargin))
    }
  }

  private def absolutizeOutputWdlFile(value: WdlValue, cwd: Path): WdlValue = value match {
    case f: WdlFile if f.valueString(0) != '/' => WdlFile(cwd.resolve(f.valueString).toAbsolutePath.toString)
    case a: WdlArray => a map {absolutizeOutputWdlFile(_, cwd)}
    case x => x
  }
}
