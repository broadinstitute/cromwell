package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Files, Path, Paths}

import better.files.{File => ScalaFile, _}
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import wdl4s.{Call, TaskOutput}
import wdl4s.types.{WdlArrayType, WdlFileType, WdlMapType}
import wdl4s.values.{WdlValue, _}
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.backend.{CallLogs, LocalFileSystemBackendCall, _}
import cromwell.engine.io.IoInterface
import cromwell.engine.io.gcs.{GcsPath, GoogleCloudStorage}
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.engine._
import cromwell.util.TryUtil
import org.apache.commons.io.FileUtils
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import Hashing._

object SharedFileSystem {
  type LocalizationStrategy = (String, Path, WorkflowDescriptor) => Try[Unit]

  val SharedFileSystemConf = ConfigFactory.load.getConfig("backend").getConfig("shared-filesystem")
  val CromwellExecutionRoot = SharedFileSystemConf.getString("root")
  val LocalizationStrategies = SharedFileSystemConf.getStringList("localization").asScala
  lazy val gcsInterface = GoogleCloudStorage.cromwellAuthenticated

  val Localizers = localizePathAlreadyLocalized _ +: (LocalizationStrategies map {
    case "hard-link" => localizePathViaHardLink _
    case "soft-link" => localizePathViaSymbolicLink _
    case "copy" => localizePathViaCopy _
    case unsupported => throw new UnsupportedOperationException(s"Localization strategy $unsupported is not recognized")
  }).+:(localizeFromGcs _)

  // Note that any unrecognized configuration will be raised when Localizers (just above) gets resolved.
  val DockerLocalizers = localizePathAlreadyLocalized _ +: (LocalizationStrategies collect {
    case "hard-link" => localizePathViaHardLink _
    case "copy" => localizePathViaCopy _
  }).+:(localizeFromGcs _)

  private def localizeFromGcs(originalPath: String, executionPath: Path, descriptor: WorkflowDescriptor): Try[Unit] = Try {
    import PathString._
    assert(originalPath.isGcsUrl)
    val content = descriptor.gcsInterface.get.downloadObject(GcsPath(originalPath))
    new ScalaFile(executionPath).createIfNotExists().write(content)
  }

  /**
    * Return a `Success` result if the file has already been localized, otherwise `Failure`.
    */
  private def localizePathAlreadyLocalized(originalPath: String, executionPath: Path, descriptor: WorkflowDescriptor): Try[Unit] = {
    if (Files.exists(executionPath)) Success(()) else Failure(new Throwable)
  }

  private def localizePathViaCopy(originalPath: String, executionPath: Path, descriptor: WorkflowDescriptor): Try[Unit] = {
    Try(Paths.get(originalPath)) map { srcPath =>
      if (Files.isDirectory(srcPath)) {
        FileUtils.copyDirectory(srcPath.toFile, executionPath.toFile)
      } else {
        FileUtils.copyFile(srcPath.toFile, executionPath.toFile)
      }
    }
  }

  private def localizePathViaHardLink(originalPath: String, executionPath: Path, descriptor: WorkflowDescriptor): Try[Unit] =
    Try(Files.createLink(executionPath, Paths.get(originalPath)))

  /**
    * TODO: The 'call' parameter here represents the call statement in WDL that references this path.
    * We're supposed to not use symbolic links if the call uses Docker.  However, this is currently a
    * bit incorrect because multiple calls can reference the same path if that path is in a declaration.
    *
    * The symbolic link will only fail in the Docker case if a Call uses the file directly and not
    * indirectly through one of its input expressions
    */

  private def localizePathViaSymbolicLink(originalPath: String, executionPath: Path, descriptor: WorkflowDescriptor): Try[Unit] = {
    Try(Paths.get(originalPath)) map { srcPath =>
      if (srcPath.toFile.isDirectory)
        Failure(new UnsupportedOperationException("Cannot localize directory with symbolic links"))
      else Files.createSymbolicLink(executionPath, srcPath.toAbsolutePath)
    }
  }

  val sharedFsFileHasher: FileHasher = { wdlFile: WdlFile => SymbolHash(ScalaFile(wdlFile.value).md5) }
}

trait SharedFileSystem {
  import SharedFileSystem._

  def useCachedCall(cachedBackendCall: LocalFileSystemBackendCall, backendCall: LocalFileSystemBackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    val source = cachedBackendCall.callRootPath.toAbsolutePath.toString
    val dest = backendCall.callRootPath.toAbsolutePath.toString
    val outputs = for {
      _ <- Try(backendCall.workflowDescriptor.ioManager.copy(source, dest))
      outputs <- postProcess(backendCall)
    } yield outputs

    outputs match {
      case Success(o) =>
        cachedBackendCall.hash map { h => CompletedExecutionHandle(SuccessfulExecution(o, Seq.empty[ExecutionEventEntry], cachedBackendCall.returnCode.contentAsString.stripLineEnd.toInt, h, Option(cachedBackendCall))) }
      case Failure(ex) => FailedExecutionHandle(ex).future
    }
  } flatten

  def workflowContext(workflowOptions: WorkflowOptions, workflowId: WorkflowId, name: String): WorkflowContext = {
    new WorkflowContext(LocalBackend.hostExecutionPath(name, workflowId).toString)
  }

  def engineFunctions(ioInterface: IoInterface, workflowContext: WorkflowContext): WorkflowEngineFunctions = {
    new LocalWorkflowEngineFunctions(ioInterface, workflowContext)
  }

  def postProcess(backendCall: LocalFileSystemBackendCall): Try[CallOutputs] = {
    implicit val hasher = backendCall.workflowDescriptor.fileHasher
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
      val unwrappedMap = outputMappings collect { case (name, Success(wdlValue)) =>
        name -> CallOutput(wdlValue, wdlValue.getHash(backendCall.workflowDescriptor))
      }
      Success(unwrappedMap.toMap)
    } else {
      val message = taskOutputFailures collect { case (name, Failure(e)) => s"$name: $e\n${ExceptionUtils.getStackTrace(e)}" }
      Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: ${message.mkString("\n")}"))
    }
  }

  def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): CallLogs = {
    val dir = LocalBackend.hostCallPath(descriptor.namespace.workflow.unqualifiedName, descriptor.id, callName, index)
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
   * NOTE: This ends up being a backdoor implementation of Backend.adjustInputPaths as both LocalBackend and SgeBackend
   *    end up with this implementation and thus use it to satisfy their contract with Backend.
   *    This is yuck-tastic and I consider this a FIXME, but not for this refactor
   */
  def adjustInputPaths(callKey: CallKey,
                       runtimeAttributes: CromwellRuntimeAttributes,
                       inputs: CallInputs,
                       workflowDescriptor: WorkflowDescriptor): CallInputs = {
    import PathString._

    val strategies = if (runtimeAttributes.docker.isDefined) DockerLocalizers else Localizers

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
    def toCallPath(path: String): Path = {
      val callDirectory = LocalBackend.hostCallPath(workflowDescriptor, callKey.scope.unqualifiedName, callKey.index)
      // Concatenate call directory with absolute input path
      val localInputPath = if(path.isGcsUrl) {
        val gcsPath = GcsPath(path)
        s"${gcsPath.bucket}/${gcsPath.objectName}"
      } else {
        Paths.get(path).toAbsolutePath.toString
      }
      Paths.get(callDirectory.toAbsolutePath.toString, localInputPath)
    }

    // Optional function to adjust the path to "docker path" if the call runs in docker
    val postProcessor: Option[Path => Path] = runtimeAttributes.docker map { _ => toDockerPath _ }
    val localizeFunction = localizeWdlValue(workflowDescriptor, toCallPath, strategies.toStream, postProcessor) _
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
  def localizeWdlValue(descriptor: WorkflowDescriptor, toDestPath: (String => Path), strategies: Stream[LocalizationStrategy], postProcessor: Option[Path => Path] = None)(wdlValue: WdlValue): Try[WdlValue] = {

    def localize(source: String, dest: Path) = strategies map { _(source, dest, descriptor) } find { _.isSuccess } getOrElse {
      Failure(new UnsupportedOperationException(s"Could not localize $source -> $dest"))
    }

    def adjustArray(t: WdlArrayType, inputArray: Seq[WdlValue]): Try[WdlArray] = {
      val tryAdjust = inputArray map localizeWdlValue(descriptor, toDestPath, strategies, postProcessor)

      TryUtil.sequence(tryAdjust, s"Failed to localize files in input Array ${wdlValue.valueString}") map { adjusted =>
        new WdlArray(t, adjusted)
      }
    }

    def adjustMap(t: WdlMapType, inputMap: Map[WdlValue, WdlValue]): Try[WdlMap] = {
      val tryAdjust = inputMap mapValues { localizeWdlValue(descriptor, toDestPath, strategies, postProcessor) }

      TryUtil.sequenceMap(tryAdjust, s"Failed to localize files in input Map ${wdlValue.valueString}") map { adjusted =>
        new WdlMap(t, adjusted)
      }
    }

    def adjustFile(path: String) = {
      val adjustedPath = toDestPath(path)
      localize(path, adjustedPath) map { Unit =>
        val finalPath = postProcessor map { _(adjustedPath) } getOrElse adjustedPath
        WdlFile(finalPath.toAbsolutePath.toString)
      }
    }

    wdlValue match {
      case wdlFile: WdlFile => adjustFile(wdlFile.value)
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
