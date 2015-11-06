package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Files, Path, Paths}

import better.files.{File => ScalaFile}
import com.typesafe.config.ConfigFactory
import cromwell.binding._
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlMapType}
import cromwell.binding.values.{WdlValue, _}
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.{ExecutionHandle, BackendCall, CallLogs, LocalFileSystemBackendCall}
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.util.TryUtil
import org.apache.commons.io.FileUtils

import scala.collection.JavaConverters._
import scala.concurrent.{Future, ExecutionContext}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystem {
  type LocalizationStrategy = (Path, Path) => Try[Unit]

  val SharedFileSystemConf = ConfigFactory.load.getConfig("backend").getConfig("shared-filesystem")
  val CromwellExecutionRoot = SharedFileSystemConf.getString("root")
  val LocalizationStrategies = SharedFileSystemConf.getStringList("localization").asScala

  val Localizers = localizePathAlreadyLocalized _ +: (LocalizationStrategies map {
    case "hard-link" => localizePathViaHardLink _
    case "soft-link" => localizePathViaSymbolicLink _
    case "copy" => localizePathViaCopy _
    case unsupported => throw new UnsupportedOperationException(s"Localization strategy $unsupported is not recognized")
  })

  // Note that any unrecognized configuration will be raised when Localizers (just above) gets resolved.
  val DockerLocalizers = localizePathAlreadyLocalized _ +: (LocalizationStrategies collect {
    case "hard-link" => localizePathViaHardLink _
    case "copy" => localizePathViaCopy _
  })

  /**
   * Return a `Success` result if the file has already been localized, otherwise `Failure`.
   */
  private def localizePathAlreadyLocalized(originalPath: Path, executionPath: Path): Try[Unit] = {
    if (Files.exists(executionPath)) Success(()) else Failure(new Throwable)
  }

  private def localizePathViaCopy(originalPath: Path, executionPath: Path): Try[Unit] = {
    if (Files.isDirectory(originalPath)) {
      Try(FileUtils.copyDirectory(originalPath.toFile, executionPath.toFile))
    } else {
      Try(FileUtils.copyFile(originalPath.toFile, executionPath.toFile))
    }
  }

  private def localizePathViaHardLink(originalPath: Path, executionPath: Path): Try[Unit] =
    Try(Files.createLink(executionPath, originalPath))

  /**
   * TODO: The 'call' parameter here represents the call statement in WDL that references this path.
   * We're supposed to not use symbolic links if the call uses Docker.  However, this is currently a
   * bit incorrect because multiple calls can reference the same path if that path is in a declaration.
   *
   * The symbolic link will only fail in the Docker case if a Call uses the file directly and not
   * indirectly through one of its input expressions
   */
  private def localizePathViaSymbolicLink(originalPath: Path, executionPath: Path): Try[Unit] = {
    if (originalPath.toFile.isDirectory)
      Failure(new UnsupportedOperationException("Cannot localize directory with symbolic links"))
    else Try(Files.createSymbolicLink(executionPath, originalPath.toAbsolutePath))
  }

  val sharedFSFileHasher: FileHasher = { wdlFile: WdlFile => SymbolHash(ScalaFile(wdlFile.value).md5) }
}

class SharedFileSystemIOInterface extends IOInterface {
  import better.files._

  override def readFile(path: String): String = Paths.get(path).contentAsString

  override def writeFile(path: String, content: String): Unit = Paths.get(path).write(content)

  override def listContents(path: String): Iterable[String] = Paths.get(path).list map { _.path.toAbsolutePath.toString } toIterable

  override def exists(path: String): Boolean = Paths.get(path).exists

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = {
    val file = Files.createTempFile(Paths.get(path), prefix, suffix)
    file.write(content)
    file.fullPath
  }

  override def glob(path: String, pattern: String): Seq[String] = {
    Paths.get(path).glob(pattern) map { _.path.fullPath } toSeq
  }

  override def copy(from: String, to: String): Unit = Paths.get(from).copyTo(Paths.get(to))
}

trait SharedFileSystem {

  import SharedFileSystem._

  type IOInterface = SharedFileSystemIOInterface

  def engineFunctions(interface: IOInterface): WdlStandardLibraryFunctions = new LocalEngineFunctionsWithoutCallContext(interface)
  def fileHasher(workflow: WorkflowDescriptor) = sharedFSFileHasher

  def ioInterface(workflowOptions: WorkflowOptions): IOInterface = new SharedFileSystemIOInterface

  def postProcess(backendCall: LocalFileSystemBackendCall): Try[CallOutputs] = {
    implicit val hasher = fileHasher(backendCall.workflowDescriptor)
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
      val unwrappedMap = outputMappings collect { case (name, Success(wdlValue)) => name -> CallOutput(wdlValue, wdlValue.getHash) }
      Success(unwrappedMap.toMap)
    } else {
      val message = taskOutputFailures collect { case (name, Failure(e)) => s"$name: $e" }
      Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: ${message.mkString("\n")}"))
    }
  }

  def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): CallLogs = {
    val dir = LocalBackend.hostCallPath(descriptor.namespace.workflow.name, descriptor.id, callName, index)
    CallLogs(
      stdout = WdlFile(dir.resolve("stdout").toAbsolutePath.toString),
      stderr = WdlFile(dir.resolve("stderr").toAbsolutePath.toString)
    )
  }
  /**
   * Creates host execution directory.
   */
  def initializeForWorkflow(descriptor: WorkflowDescriptor): Try[HostInputs] = {
    val hostExecutionDirectory = LocalBackend.hostExecutionPath(descriptor).toFile
    hostExecutionDirectory.mkdirs()
    Success(descriptor.actualInputs)
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(callKey: CallKey, inputs: CallInputs, workflowDescriptor: WorkflowDescriptor): CallInputs = {
    val call = callKey.scope
    val strategies = if (call.docker.isDefined) DockerLocalizers else Localizers

    def toDockerPath(path: Path): Path = {
      // Host path would look like cromwell-executions/three-step/f00ba4/call-ps/stdout.txt
      // Container path should look like /root/f00ba4/call-ps/stdout.txt
      val fullPath = path.toFile.getAbsolutePath
      // Strip out everything before cromwell-executions.
      val pathUnderCromwellExecutions = fullPath.substring(fullPath.indexOf(CromwellExecutionRoot) + CromwellExecutionRoot.length)
      // Strip out the workflow name (the first component under cromwell-executions).
      val pathWithWorkflowName = Paths.get(pathUnderCromwellExecutions)
      Paths.get(WdlFile.appendPathsWithSlashSeparators("/root", pathWithWorkflowName.subpath(1, pathWithWorkflowName.getNameCount).toString))
    }

    /**
     * Transform an original input path to a path in the call directory.
     * The new path matches the original path, it only "moves" the root to be the call directory.
     */
    def toCallPath(path: Path): Path = {
      val callDirectory = LocalBackend.hostCallPath(workflowDescriptor, call.name, callKey.index)
      // Concatenate call directory with absolute input path
      Paths.get(callDirectory.toAbsolutePath.toString, path.toAbsolutePath.toString)
    }

    // Optional function to adjust the path to "docker path" if the call runs in docker
    val postProcessor: Option[Path => Path] = call.docker map { _ => toDockerPath _ }
    val localizeFunction = localizeWdlValue(toCallPath, strategies.toStream, postProcessor) _
    val localizedValues = inputs.toSeq map {
      case (name, value) => localizeFunction(value) map { name -> _ }
    }

    TryUtil.sequence(localizedValues, "Failures during localization").get toMap
  }

  /**
   * Try to localize a WdlValue if it is or contains a WdlFile.
   * @param toDestPath function specifying how to generate the destination path from the source path
   * @param strategies strategies to use for localization
   * @param postProcessor optional function to be applied to the path after the file it points to has been localized (defaults to noop)
   * @param wdlValue WdlValue to localize
   * @return localized wdlValue
   */
  def localizeWdlValue(toDestPath: (Path => Path), strategies: Stream[LocalizationStrategy], postProcessor: Option[Path => Path] = None)(wdlValue: WdlValue): Try[WdlValue] = {

    def localize(source: Path, dest: Path) = strategies map { _(source, dest) } find { _.isSuccess } getOrElse {
      Failure(new UnsupportedOperationException(s"Could not localize $source -> $dest"))
    }

    def adjustArray(t: WdlArrayType, inputArray: Seq[WdlValue]): Try[WdlArray] = {
      val tryAdjust = inputArray map localizeWdlValue(toDestPath, strategies, postProcessor)

      TryUtil.sequence(tryAdjust, s"Failed to localize files in input Array ${wdlValue.valueString}") map { adjusted =>
        new WdlArray(t, adjusted)
      }
    }

    def adjustMap(t: WdlMapType, inputMap: Map[WdlValue, WdlValue]): Try[WdlMap] = {
      val tryAdjust = inputMap mapValues { localizeWdlValue(toDestPath, strategies, postProcessor) }

      TryUtil.sequenceMap(tryAdjust, s"Failed to localize files in input Map ${wdlValue.valueString}") map { adjusted =>
        new WdlMap(t, adjusted)
      }
    }

    def adjustFile(path: Path) = {
      val adjustedPath = toDestPath(path)
      localize(path, adjustedPath) map { Unit =>
        val finalPath = postProcessor map { _(adjustedPath) } getOrElse adjustedPath
        WdlFile(finalPath.toAbsolutePath.toString)
      }
    }

    wdlValue match {
      case wdlFile: WdlFile => adjustFile(Paths.get(wdlFile.value))
      case WdlArray(t, values) => adjustArray(t, values)
      case WdlMap(t, values) => adjustMap(t, values)
      case x => Success(x)
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
      case rhs if rhs.wdlType == taskOutput.wdlType => Success(rhs)
      case rhs: WdlString if taskOutput.wdlType == WdlFileType => assertTaskOutputPathExists(hostAbsoluteFilePath(backendCall, rhs.value), taskOutput, backendCall.call.fullyQualifiedName)
      case rhs => taskOutput.wdlType.coerceRawValue(rhs)
    }
  }

  private def absolutizeOutputWdlFile(value: WdlValue, cwd: Path): WdlValue = value match {
    case f: WdlFile if f.valueString(0) != '/' => WdlFile(cwd.resolve(f.valueString).toAbsolutePath.toString)
    case a: WdlArray => a map {absolutizeOutputWdlFile(_, cwd)}
    case x => x
  }
}
