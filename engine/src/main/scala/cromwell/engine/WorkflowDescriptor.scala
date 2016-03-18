package cromwell.engine

import java.nio.file._

import better.files._
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.engine.backend.io._
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.backend.{Backend, BackendType, CromwellBackend}
import cromwell.engine.workflow.WorkflowOptions
import cromwell.logging.WorkflowLogger
import cromwell.util.{SimpleExponentialBackoff, TryUtil}
import cromwell.webservice.WorkflowMetadataResponse
import lenthall.config.ScalaConfig._
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}
import spray.json.{JsObject, _}
import wdl4s._
import wdl4s.values.{WdlFile, WdlSingleFile, WdlValue, _}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._

case class WorkflowDescriptor(id: WorkflowId,
                              sourceFiles: WorkflowSourceFiles,
                              workflowOptions: WorkflowOptions,
                              workflowLogOptions: Option[WorkflowLogOptions],
                              rawInputs: Map[String, JsValue],
                              namespace: NamespaceWithWorkflow,
                              coercedInputs: WorkflowCoercedInputs,
                              declarations: WorkflowCoercedInputs,
                              backend: Backend,
                              configCallCaching: Boolean,
                              lookupDockerHash: Boolean,
                              wfContext: WorkflowContext,
                              fileSystems: List[FileSystem]) {
  import WorkflowDescriptor._

  val shortId = id.toString.split("-")(0)
  val name = namespace.workflow.unqualifiedName
  val actualInputs: WorkflowCoercedInputs = coercedInputs ++ declarations
  private val relativeWorkflowRootPath = s"$name/$id"
  private val log = WorkflowLogger("WorkflowDescriptor", this)
  val workflowLogDir = workflowOptions.get(WorkflowLogDirOptionKey) recover { case e: IllegalArgumentException =>
    log.warn(s"$WorkflowLogDirOptionKey expected to be of type String", e)
    throw e
  }
  val workflowOutputsPath = workflowOptions.get(WorkflowOutputsOptionKey) recover { case e: IllegalArgumentException =>
    log.warn(s"$WorkflowOutputsOptionKey expected to be of type String", e)
    throw e
  }
  val callLogsDir = workflowOptions.get(CallLogsDirOptionKey) recover { case e: IllegalArgumentException =>
    log.warn(s"$CallLogsDirOptionKey expected to be of type String", e)
    throw e
  }

  lazy val fileHasher: FileHasher = { wdlFile: WdlFile =>
    try {
      SymbolHash(wdlFile.value.toPath(fileSystems).hash)
    } catch {
      case e: Exception =>
        log.error(s"Cannot compute hash for file $wdlFile")
        throw e
    }
  }

  private lazy val optionCacheWriting = workflowOptions getBoolean "write_to_cache" getOrElse configCallCaching
  private lazy val optionCacheReading = workflowOptions getBoolean "read_from_cache" getOrElse configCallCaching

  if (!configCallCaching) {
    if (optionCacheWriting) logWriteDisabled()
    if (optionCacheReading) logReadDisabled()
  }

  lazy val writeToCache = configCallCaching && optionCacheWriting
  lazy val readFromCache = configCallCaching && optionCacheReading

  lazy val workflowLogName = s"workflow.$id.log"
  lazy val workflowLogPath = workflowLogOptions.map(_.dir.createDirectories() / workflowLogName).map(_.path)
  lazy val workflowLogger = workflowLogPath match {
    case Some(path) => makeFileLogger(path, workflowLogName, Level.toLevel(sys.props.getOrElse("LOG_LEVEL", "debug")))
    case None => NOPLogger.NOP_LOGGER
  }

  def maybeDeleteWorkflowLog(): Unit = {
    for {
      opt <- workflowLogOptions if opt.temporary
      log <- workflowLogPath
    } yield log.delete(ignoreIOExceptions = true)
  }


  lazy val workflowRootPath = wfContext.root.toPath(fileSystems)
  def workflowRootPathWithBaseRoot(rootPath: String): Path = {
    backend.buildWorkflowRootPath(rootPath, name, id).toPath(fileSystems)
  }

  def copyWorkflowOutputs(workflowMetadataResponse: WorkflowMetadataResponse)
                         (implicit executionContext: ExecutionContext): Future[Unit] = {
    // Try to copy outputs to final destination
    workflowOutputsPath map copyWorkflowOutputsTo(workflowMetadataResponse) getOrElse Future.successful(())
  }

  private def copyWorkflowOutputsTo(workflowMetadataResponse: WorkflowMetadataResponse)(destDirectory: String)
                             (implicit executionContext: ExecutionContext): Future[Unit] = {
    Future(copyWdlFilesTo(destDirectory, workflowMetadataResponse.outputs.toSeq.flatMap(_.values)))
  }

  def copyCallLogs(workflowMetadataResponse: WorkflowMetadataResponse)
                  (implicit executionContext: ExecutionContext): Future[Unit] = {
    callLogsDir map copyCallLogsTo(workflowMetadataResponse) getOrElse Future.successful(())
  }

  private def copyCallLogsTo(workflowMetadataResponse: WorkflowMetadataResponse)(destDirectory: String)
                            (implicit executionContext: ExecutionContext): Future[Unit] = {
    Future {
      val callLogs = for {
        callMetadatas <- workflowMetadataResponse.calls.values
        callMetadata <- callMetadatas
        callLog <- callMetadata.stdout.toSeq ++ callMetadata.stderr.toSeq ++
          callMetadata.backendLogs.toSeq.flatMap(_.values)
      } yield callLog
      copyWdlFilesTo(destDirectory, callLogs)
    }
  }

  def copyWorkflowLog()(implicit executionContext: ExecutionContext): Future[Unit] = {
    (workflowLogDir, workflowLogPath) match {
      case (Success(dir), Some(path)) =>
        val logger = backend.workflowLogger(this)
        val dest = dir.toPath(fileSystems).resolve(workflowLogName)
        Future.fromTry(copyFile(logger, dest, path))
      case _ => Future.successful(())
    }
  }

  private def copyWdlFilesTo(destDirectory: String, wdlValues: Traversable[WdlValue]): Unit = {
    val logger = backend.workflowLogger(this)

    // All outputs should have wdl values at this point, if they don't there's nothing we can do here
    val copies = for {
      wdlValue <- wdlValues
      wdlFile <- wdlValue collectAsSeq { case f: WdlSingleFile => f }
    } yield copyWdlFile(logger, destDirectory, wdlFile)

    // Throw an exception if one or more of the copies failed.
    TryUtil.sequence(copies.toSeq) match {
      case Success(_) => ()
      case Failure(e) => throw new Exception(s"Copy failed for the following files:\n ${e.getMessage}", e)
    }
  }

  def copyWdlFile(logger: WorkflowLogger, destDirectory: String, file: WdlFile): Try[Unit] = {

    val src = file.valueString.toPath(fileSystems)
    val wfPath = wfContext.root.toPath(fileSystems).toAbsolutePath
    val srcSubPath = src.subpath(wfPath.getNameCount, src.getNameCount)
    val dest = destDirectory
      .toPath(fileSystems)
      .toAbsolutePath
      .resolve(relativeWorkflowRootPath)
      // UnixPath.resolve(NioGcsPath) seems to be returning a null pointer. TODO: Add a test to confirm
      .resolve(srcSubPath.toString)
    copyFile(logger, dest, src)
  }

  def copyFile(logger: WorkflowLogger, dest: Path, src: Path): Try[Unit] = {
    def copy(): Unit = {
      logger.info(s"Trying to copy output file $src to $dest")
      Files.createDirectories(dest.getParent)
      Files.copy(src, dest)
    }

    TryUtil.retryBlock(
      fn = (retries: Option[Unit]) => copy(),
      retryLimit = Option(5),
      backoff = SimpleExponentialBackoff(5 seconds, 10 seconds, 1.1D),
      logger = logger,
      failMessage = Option(s"Failed to copy file $src to $dest"),
      isFatal = (t: Throwable) => t.isInstanceOf[FileAlreadyExistsException]
    ) recover {
      case _: FileAlreadyExistsException =>
        logger.info(s"Tried to copy the same file multiple times. Skipping subsequent copies for $src")
    }
  }

  private def makeFileLogger(path: Path, name: String, level: Level): Logger = {
    val ctx = LoggerFactory.getILoggerFactory.asInstanceOf[LoggerContext]

    /*
    WorkflowDescriptor.copy() is currently invoked by WorkflowActor.startActor().
    This causes this block of code to be executed twice.
     */
    Option(ctx.exists(name)) match {
      case Some(existingLogger) => existingLogger
      case None =>
        val encoder = new PatternLayoutEncoder()
        encoder.setPattern("%date %-5level - %msg%n")
        encoder.setContext(ctx)
        encoder.start()

        val appender = new FileAppender[ILoggingEvent]()
        appender.setFile(path.fullPath)
        appender.setEncoder(encoder)
        appender.setName(name)
        appender.setContext(ctx)
        appender.start()

        val fileLogger = ctx.getLogger(name)
        fileLogger.addAppender(appender)
        fileLogger.setAdditive(false)
        fileLogger.setLevel(level)
        fileLogger
    }
  }

  private def logWriteDisabled() = workflowLogger.warn(writeDisabled)
  private def logReadDisabled() = workflowLogger.warn(readDisabled)

  override def toString = s"WorkflowDescriptor(${id.id.toString})"

  def hash(wdlValue: WdlValue): Option[SymbolHash] = {
    if (configCallCaching) Option(wdlValue.computeHash(fileHasher)) else None
  }
}

object WorkflowDescriptor {

  val WorkflowLogDirOptionKey = "workflow_log_dir"
  val WorkflowOutputsOptionKey = "outputs_path"
  val CallLogsDirOptionKey = "call_logs_dir"
  val OptionKeys: Set[String] = Set(WorkflowLogDirOptionKey, WorkflowOutputsOptionKey, CallLogsDirOptionKey)

  def apply(id: WorkflowId, sourceFiles: WorkflowSourceFiles): WorkflowDescriptor = {
    WorkflowDescriptor(id, sourceFiles, ConfigFactory.load)
  }

  def apply(id: WorkflowId, sourceFiles: WorkflowSourceFiles, conf: Config): WorkflowDescriptor = {
    validateWorkflowDescriptor(id, sourceFiles, CromwellBackend.backend(), conf) match {
      case scalaz.Success(w) => w
      case scalaz.Failure(f) =>
        throw new IllegalArgumentException() with ThrowableWithErrors {
          val message = s"Workflow input processing failed."
          val errors = f
        }
    }
  }

  private def validateWorkflowDescriptor(id: WorkflowId,
                                         sourceFiles: WorkflowSourceFiles,
                                         backend: Backend,
                                         conf: Config): ErrorOr[WorkflowDescriptor] = {
    val namespace = validateNamespace(id, sourceFiles.wdlSource)
    val options = validateWorkflowOptions(id, sourceFiles.workflowOptionsJson, backend)

    (namespace |@| options) { (_, _) } flatMap { case (nam, opt) =>
      val runtimeAttributes = validateRuntimeAttributes(id, nam, backend.backendType)
      val rawInputs = validateRawInputs(id, sourceFiles.inputsJson)
      (runtimeAttributes |@| rawInputs) { (_, _) } flatMap { case (_, raw) =>
        buildWorkflowDescriptor(id, sourceFiles, nam, raw, opt, backend, conf)
      }
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFiles,
                                      namespace: NamespaceWithWorkflow,
                                      rawInputs: Map[String, JsValue],
                                      options: WorkflowOptions,
                                      backend: Backend,
                                      conf: Config): ErrorOr[WorkflowDescriptor] = {
    val workflowRootPath = backend.buildWorkflowRootPath(backend.rootPath(options), namespace.workflow.unqualifiedName, id)
    val wfContext = new WorkflowContext(workflowRootPath)
    val fileSystems = backend.fileSystems(options)
    val engineFunctions = backend.engineFunctions(fileSystems, wfContext)

    val validatedDescriptor = for {
      c <- validateCoercedInputs(id, rawInputs, namespace).disjunction
      d <- validateDeclarations(id, namespace, options, c, engineFunctions).disjunction
    } yield WorkflowDescriptor(id, sourceFiles, options, workflowLogOptions(conf), rawInputs, namespace, c, d, backend, configCallCaching(conf), lookupDockerHash(conf),
      wfContext, fileSystems)

    validatedDescriptor.validation
  }

  private def validateNamespace(id: WorkflowId, source: WdlSource): ErrorOr[NamespaceWithWorkflow] = {
    try {
      NamespaceWithWorkflow.load(source).successNel
    } catch {
      case e: Exception => s"Workflow $id unable to load namespace: ${e.getMessage}".failureNel
    }
  }

  private def validateRuntimeAttributes(id: WorkflowId, namespace: NamespaceWithWorkflow, backendType: BackendType): ErrorOr[Unit] = {
    Try(namespace.workflow.calls.map(_.task.runtimeAttributes) foreach { r => CromwellRuntimeAttributes.validateKeys(r.attrs.keySet, backendType) }) match {
      case scala.util.Success(_) => ().successNel
      case scala.util.Failure(e) => s"Workflow $id contains bad runtime attributes: ${e.getMessage}".failureNel
    }
  }

  private def validateWorkflowOptions(id: WorkflowId,
                                      optionsJson: WorkflowOptionsJson,
                                      backend: Backend): ErrorOr[WorkflowOptions] = {
    WorkflowOptions.fromJsonString(optionsJson) match {
      case Success(o) => validateBackendOptions(id, o, backend)
      case Failure(e) => s"Workflow $id contains bad options JSON: ${e.getMessage}".failureNel
    }
  }

  private def validateRawInputs(id: WorkflowId, json: WdlJson): ErrorOr[Map[String, JsValue]] = {
    Try(json.parseJson) match {
      case Success(JsObject(inputs)) => inputs.successNel
      case _ => s"Workflow $id contains bad inputs JSON: $json".failureNel
    }
  }

  private def validateCoercedInputs(id: WorkflowId,
                                    rawInputs: Map[String, JsValue],
                                    namespace: NamespaceWithWorkflow): ErrorOr[WorkflowCoercedInputs] = {
    namespace.coerceRawInputs(rawInputs) match {
      case Success(r) => r.successNel
      case Failure(e: ThrowableWithErrors) => scalaz.Failure(e.errors)
      case Failure(e) => e.getMessage.failureNel
    }
  }

  private def validateBackendOptions(id: WorkflowId, options: WorkflowOptions, backend: Backend): ErrorOr[WorkflowOptions] = {
    try {
      backend.assertWorkflowOptions(options)
      options.successNel
    } catch {
      case e: Exception =>
        s"Workflow $id has invalid options for backend ${backend.backendType}: ${e.getMessage}".failureNel
    }
  }

  private def validateDeclarations(id: WorkflowId,
                                   namespace: NamespaceWithWorkflow,
                                   options: WorkflowOptions,
                                   coercedInputs: WorkflowCoercedInputs,
                                   engineFunctions: WorkflowEngineFunctions): ErrorOr[WorkflowCoercedInputs] = {
    namespace.staticDeclarationsRecursive(coercedInputs, engineFunctions) match {
      case Success(d) => d.successNel
      case Failure(e) => s"Workflow $id has invalid declarations: ${e.getMessage}".failureNel
    }
  }

  private val DefaultCallCachingValue = false
  private val DefaultLookupDockerHash = false

  private def disabledMessage(readWrite: String, consequence: String) =
    s"""$readWrite is enabled in the workflow options but Call Caching is disabled in this Cromwell instance.
       |As a result the calls in this workflow $consequence
       """.stripMargin

  private val writeDisabled = disabledMessage("Write to Cache", "WILL NOT be cached")
  private val readDisabled = disabledMessage("Read from Cache", "WILL ALL be executed")

  private def configCallCaching(conf: Config) = lookupBooleanWithDefault(conf, "call-caching", "enabled", DefaultCallCachingValue)
  private def lookupDockerHash(conf: Config) = lookupBooleanWithDefault(conf, "call-caching", "lookup-docker-hash", DefaultLookupDockerHash)

  private def lookupBooleanWithDefault(conf: Config, stanza: String, key: String, default: Boolean) = {
    (for {
      config <- conf.getConfigOption(stanza)
      value <- config.getBooleanOption(key)
    } yield value) getOrElse default
  }

  private def workflowLogOptions(conf: Config): Option[WorkflowLogOptions] = {
    for {
      workflowConfig <- conf.getConfigOption("workflow-options")
      dir <- workflowConfig.getStringOption("workflow-log-dir") if !dir.isEmpty
      temporary <- workflowConfig.getBooleanOption("workflow-log-temporary") orElse Option(true)
    } yield WorkflowLogOptions(Paths.get(dir), temporary)
  }
}

case class WorkflowLogOptions(dir: Path, temporary: Boolean)
