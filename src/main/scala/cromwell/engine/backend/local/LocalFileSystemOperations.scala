package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Path, Files, Paths}

import cromwell.binding._
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlMapType}
import cromwell.binding.values.{WdlValue, _}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.WorkflowId
import cromwell.engine.backend.{StdoutStderr, LocalFileSystemBackendCall}
import org.apache.commons.io.FileUtils

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

trait LocalFileSystemOperations {
  def postProcess(backendCall: LocalFileSystemBackendCall): Try[CallOutputs] = {
    // Evaluate output expressions, performing conversions from String -> File where required.
    val outputMappings = backendCall.call.task.outputs map { taskOutput =>
      val tryConvertedValue =
        for {
          expressionValue <- taskOutput.expression.evaluate(backendCall.lookupFunction, backendCall.engineFunctions, interpolateStrings=true)
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

  def stdoutStderr(workflowId: WorkflowId, workflowName: String, callName: String, index: ExecutionIndex): StdoutStderr = {
    val dir = LocalBackend.hostCallPath(workflowName, workflowId, callName, index)
    StdoutStderr(
      stdout = WdlFile(dir.resolve("stdout").toAbsolutePath.toString),
      stderr = WdlFile(dir.resolve("stderr").toAbsolutePath.toString)
    )
  }

  /**
   * Creates host execution directory, inputs path, and outputs path.  Stages any input files into the workflow-inputs
   * directory and localizes their paths relative to the container.
   */
  def initializeForWorkflow(descriptor: WorkflowDescriptor): HostInputs = {
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
      val pathUnderCromwellExecutions = fullPath.substring(fullPath.indexOf(LocalBackend.CromwellExecutions) + LocalBackend.CromwellExecutions.length)
      // Strip out the workflow name (the first component under cromwell-executions).
      val pathWithWorkflowName = Paths.get(pathUnderCromwellExecutions)
      WdlFile(WdlFile.appendPathsWithSlashSeparators("/root", pathWithWorkflowName.subpath(1, pathWithWorkflowName.getNameCount).toString))
    }

    // If this call is using Docker, adjust input paths, otherwise return unaltered input paths.
    def adjustPath(nameAndValue: (String, WdlValue)): (String, WdlValue) = {
      val (name, value) = nameAndValue
      val adjusted = value match {
        case WdlFile(path) => containerPath(path)
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
  private def stageWorkflowInputs(descriptor: WorkflowDescriptor): HostInputs = {
    val hostInputsPath = Paths.get(LocalBackend.hostExecutionPath(descriptor).toFile.getAbsolutePath, "workflow-inputs")
    descriptor.actualInputs map {case(name, value) => name -> stageWdlValue(value, hostInputsPath)}
  }

  private def stageInputArray(array: WdlArray, hostInputsPath: Path): WdlArray = array.map(stageWdlValue(_, hostInputsPath))

  private def stageWdlValue(value: WdlValue, hostInputsPath: Path): WdlValue = value match {
    case w: WdlFile => stageWdlFile(w, hostInputsPath)
    case a: WdlArray => stageInputArray(a, hostInputsPath)
    case x => x
  }

  private def stageWdlFile(wdlFile: WdlFile, hostInputsPath: Path): WdlFile = {
    val originalPath = Paths.get(wdlFile.value)
    val executionPath = hostInputsPath.resolve(originalPath.getFileName.toString)
    if (Files.isDirectory(originalPath)) {
      FileUtils.copyDirectory(originalPath.toFile, executionPath.toFile)
    } else {
      FileUtils.copyFile(originalPath.toFile, executionPath.toFile)
    }
    WdlFile(executionPath.toString)
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
