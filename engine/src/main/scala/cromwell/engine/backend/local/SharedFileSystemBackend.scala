package cromwell.engine.backend.local

import java.nio.file._

import better.files._
import cromwell.WorkflowEngineFunctions
import cromwell.core.{CallOutput, CallOutputs, WorkflowOptions, _}
import cromwell.engine.backend._
import cromwell.engine.io.gcs.GcsPath
import cromwell.engine.{ExecutionEventEntry, backend}
import cromwell.util.TryUtil
import org.apache.commons.lang3.exception.ExceptionUtils
import wdl4s.types.{WdlArrayType, WdlFileType, WdlMapType}
import wdl4s.values.{WdlValue, _}
import wdl4s.{Call, CallInputs, TaskOutput}

import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object SharedFileSystemBackend {
  type LocalizationStrategy = (String, Path, WorkflowDescriptor) => Try[Unit]

  private[local] def localizeFromGcs(originalPath: String, executionPath: Path,
                                     descriptor: WorkflowDescriptor): Try[Unit] = Try {
    import backend.io._
    assert(originalPath.isGcsUrl)
    Option(executionPath.parent) map { _.createDirectories }
    originalPath.toAbsolutePath(descriptor.fileSystems).copyTo(executionPath)
  }

  /**
    * Return a `Success` result if the file has already been localized, otherwise `Failure`.
    */
  private[local] def localizePathAlreadyLocalized(originalPath: String, executionPath: Path,
                                                  descriptor: WorkflowDescriptor): Try[Unit] = {
    Try {
      if (!Files.exists(executionPath))
        throw new NoSuchFileException(s"$executionPath does not exist")
    }
  }

  private[local] def localizePathViaCopy(originalPath: String, executionPath: Path,
                                         descriptor: WorkflowDescriptor): Try[Unit] = {
    Try {
      val srcPath = File(originalPath)
      Option(executionPath.parent) map { _.createDirectories }
      srcPath.copyTo(executionPath)
    }
  }

  private[local] def localizePathViaHardLink(originalPath: String, executionPath: Path,
                                             descriptor: WorkflowDescriptor): Try[Unit] =
    Try {
      val srcPath = File(originalPath)
      Option(executionPath.parent) map { _.createDirectories }
      executionPath.linkTo(srcPath, symbolic = false)
    }

  /**
    * TODO: The 'call' parameter here represents the call statement in WDL that references this path.
    * We're supposed to not use symbolic links if the call uses Docker.  However, this is currently a
    * bit incorrect because multiple calls can reference the same path if that path is in a declaration.
    *
    * The symbolic link will only fail in the Docker case if a Call uses the file directly and not
    * indirectly through one of its input expressions
    */

  private[local] def localizePathViaSymbolicLink(originalPath: String, executionPath: Path,
                                                 descriptor: WorkflowDescriptor): Try[Unit] = {
    Try {
      val srcPath = File(originalPath)
      if (srcPath.isDirectory)
        throw new UnsupportedOperationException("Cannot localize directory with symbolic links")
      Option(executionPath.parent) map { _.createDirectories }
      executionPath.linkTo(srcPath, symbolic = true)
    }
  }

  val sharedFsFileHasher: FileHasher = { wdlFile: WdlFile => SymbolHash(File(wdlFile.value).md5) }
}

trait SharedFileSystemBackend extends CanUseGcsFilesystem { self: Backend =>
  import SharedFileSystemBackend._
  import backend.io._

  val CromwellExecutionRoot = backendConfig.getString("root")
  val LocalizationStrategies = backendConfig.getConfig("filesystems.local").getStringList("localization").asScala

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

  def useCachedCall(cachedJobDescriptor: BackendCallJobDescriptor, jobDescriptor: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    val source = cachedJobDescriptor.callRootPath
    // Use the same filesystems for the destination path as the source path so they share the same providers.
    val dest = jobDescriptor.callRootPath.toAbsolutePath.toString.toAbsolutePath(cachedJobDescriptor.workflowDescriptor.fileSystems)
    val outputs = for {
      _ <- Try(source.copyTo(dest))
      outputs <- postProcess(jobDescriptor)
    } yield outputs

    outputs match {
      case Success(o) =>
        cachedJobDescriptor.hash map { h => CompletedExecutionHandle(SuccessfulBackendCallExecution(o, Seq.empty[ExecutionEventEntry], callRootPath(jobDescriptor).resolve("rc").contentAsString.stripLineEnd.toInt, h, Option(cachedJobDescriptor))) }
      case Failure(ex) => FailedExecutionHandle(ex).future
    }
  } flatten

  def rootPath(workflowOptions: WorkflowOptions) = CromwellExecutionRoot

  def engineFunctions(fileSystems: List[FileSystem], workflowContext: WorkflowContext): WorkflowEngineFunctions = {
    new LocalWorkflowEngineFunctions(fileSystems, workflowContext)
  }

  def postProcess(jobDescriptor: BackendCallJobDescriptor): Try[CallOutputs] = {
    implicit val hasher = jobDescriptor.workflowDescriptor.fileHasher

    val outputs = jobDescriptor.call.task.outputs
    val outputFoldingFunction = getOutputFoldingFunction(jobDescriptor)
    val outputMappings = outputs.foldLeft(Seq.empty[AttemptedLookupResult])(outputFoldingFunction).map(_.toPair).toMap

    val taskOutputFailures = outputMappings filter { _._2.isFailure }
    if (taskOutputFailures.isEmpty) {
      val unwrappedMap = outputMappings collect { case (name, Success(wdlValue)) =>
        name -> CallOutput(wdlValue, jobDescriptor.workflowDescriptor.hash(wdlValue))
      }
      Success(unwrappedMap)
    } else {
      val message = taskOutputFailures collect { case (name, Failure(e)) => s"$name: $e\n${ExceptionUtils.getStackTrace(e)}" }
      Failure(new Throwable(s"Workflow ${jobDescriptor.workflowDescriptor.id}: ${message.mkString("\n")}"))
    }
  }

  private def getOutputFoldingFunction(jobDescriptor: BackendCallJobDescriptor): (Seq[AttemptedLookupResult], TaskOutput) => Seq[AttemptedLookupResult] = {
    (currentList: Seq[AttemptedLookupResult], taskOutput: TaskOutput) => {
      currentList ++ Seq(AttemptedLookupResult(taskOutput.name, outputLookup(taskOutput, jobDescriptor, currentList)))
    }
  }

  private def outputLookup(taskOutput: TaskOutput, jobDescriptor: BackendCallJobDescriptor, currentList: Seq[AttemptedLookupResult]) = for {
    expressionValue <- taskOutput.requiredExpression.evaluate(jobDescriptor.lookupFunction(currentList.toLookupMap), jobDescriptor.callEngineFunctions)
    convertedValue <- outputAutoConversion(jobDescriptor, taskOutput, expressionValue)
    pathAdjustedValue <- Success(absolutizeOutputWdlFile(convertedValue, jobDescriptor.callRootPath))
  } yield pathAdjustedValue

  def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  def sharedFileSystemStdoutStderr(jobDescriptor: BackendCallJobDescriptor): CallLogs = {
    val dir = jobDescriptor.callRootPath
    CallLogs(
      stdout = WdlFile(dir.resolve("stdout").toAbsolutePath.toString),
      stderr = WdlFile(dir.resolve("stderr").toAbsolutePath.toString)
    )
  }

  /**
   * Creates host execution directory.
   */
  def initializeForWorkflow(descriptor: WorkflowDescriptor): Try[Unit] = Try {
    val hostExecutionDirectory = descriptor.workflowRootPath.toFile
    hostExecutionDirectory.mkdirs()
  }

  def toDockerPath(path: WdlValue): WdlValue = path match {
    case file: WdlFile => WdlFile(toDockerPath(Paths.get(path.valueString)).toAbsolutePath.toString)
    case v => v
  }

  private def toDockerPath(path: Path): Path = {
    path.toAbsolutePath match {
      case p if p.startsWith(LocalBackend.ContainerRoot) => p
      case p =>
        /** For example:
          *
          * p = /abs/path/to/cromwell-executions/three-step/f00ba4/call-ps/stdout.txt
          * localExecutionRoot = /abs/path/to/cromwell-executions
          * subpath = three-step/f00ba4/call-ps/stdout.txt
          *
          * return value = /root/three-step/f00ba4/call-ps/stdout.txt
          *
          * TODO: this assumes that p.startsWith(localExecutionRoot)
          */
        val localExecutionRoot = Paths.get(CromwellExecutionRoot).toAbsolutePath
        val subpath = p.subpath(localExecutionRoot.getNameCount, p.getNameCount)
        Paths.get(LocalBackend.ContainerRoot).resolve(subpath)
    }
  }

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   * NOTE: This ends up being a backdoor implementation of Backend.adjustInputPaths as both LocalBackend and SgeBackend
   *    end up with this implementation and thus use it to satisfy their contract with Backend.
   *    This is yuck-tastic and I consider this a FIXME, but not for this refactor
   */
  def adjustSharedInputPaths(jobDescriptor: BackendCallJobDescriptor): CallInputs = {

    val strategies = if (jobDescriptor.callRuntimeAttributes.docker.isDefined) DockerLocalizers else Localizers

    /**
      * Transform an original input path to a path in the call directory.
      * The new path matches the original path, it only "moves" the root to be the call directory.
      */
    def toCallPath(path: String): Path = {
      val callDirectory = jobDescriptor.callRootPath
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
    val postProcessor: Option[Path => Path] = jobDescriptor.callRuntimeAttributes.docker map { _ => toDockerPath _ }
    val localizeFunction = localizeWdlValue(jobDescriptor.workflowDescriptor, toCallPath, strategies.toStream, postProcessor) _
    val localizedValues = jobDescriptor.locallyQualifiedInputs.toSeq map {
      case (name, value) => localizeFunction(value) map { name -> _ }
    }

    TryUtil.sequence(localizedValues, "Failures during localization").get toMap
  }

  /**
   * Try to localize a WdlValue if it is or contains a WdlFile.
   *
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

  private def hostAbsoluteFilePath(jobDescriptor: BackendCallJobDescriptor, pathString: String): String =
    if (Paths.get(pathString).isAbsolute) {
      pathString
    } else {
      Paths.get(jobDescriptor.callRootPath.fullPath, pathString).fullPath
    }

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
  private def outputAutoConversion(jobDescriptor: BackendCallJobDescriptor, taskOutput: TaskOutput, rawOutputValue: WdlValue): Try[WdlValue] = {
    rawOutputValue match {
      case rhs if rhs.wdlType == taskOutput.wdlType => Success(rhs)
      case rhs: WdlString if taskOutput.wdlType == WdlFileType => assertTaskOutputPathExists(hostAbsoluteFilePath(jobDescriptor, rhs.value), taskOutput, jobDescriptor.call.fullyQualifiedName)
      case rhs => taskOutput.wdlType.coerceRawValue(rhs)
    }
  }

  private def absolutizeOutputWdlFile(value: WdlValue, cwd: Path): WdlValue = value match {
    case f: WdlFile if f.valueString(0) != '/' => WdlFile(cwd.resolve(f.valueString).toAbsolutePath.toString)
    case a: WdlArray => a map {absolutizeOutputWdlFile(_, cwd)}
    case x => x
  }

  protected def buildCallContext(descriptor: BackendCallJobDescriptor): CallContext = {
    val callRoot = callRootPath(descriptor)
    val stdout = callRoot.resolve("stdout")
    val stderr = callRoot.resolve("stderr")

    new CallContext(callRoot.fullPath, stdout.fullPath, stderr.fullPath)
  }

  override def fileSystems(options: WorkflowOptions): List[FileSystem] = {
    // The default local filesystem should already have been validated, check the GCS filesystem.
    List(gcsFilesystem(options), Option(defaultFileSystem)).flatten
  }
}
