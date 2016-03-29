package cromwell.engine.backend

import java.nio.file._

import better.files._
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import cromwell.core.{WorkflowId, WorkflowOptions}
import cromwell.engine.{WorkflowFailureMode, WorkflowSourceFiles}
import cromwell.engine.backend.io._
import cromwell.logging.WorkflowLogger
import cromwell.util.{SimpleExponentialBackoff, TryUtil}
import cromwell.webservice.WorkflowMetadataResponse
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}
import spray.json._
import wdl4s._
import wdl4s.values.{WdlFile, WdlSingleFile, WdlValue, _}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class WorkflowLogOptions(dir: Path, temporary: Boolean)

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
                              workflowFailureMode: WorkflowFailureMode,
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
  val WorkflowFailureModeKey = "workflowFailureMode"
  val OptionKeys: Set[String] = Set(WorkflowLogDirOptionKey, WorkflowOutputsOptionKey, CallLogsDirOptionKey, WorkflowFailureModeKey)

  private def disabledMessage(readWrite: String, consequence: String) =
    s"""$readWrite is enabled in the workflow options but Call Caching is disabled in this Cromwell instance.
       |As a result the calls in this workflow $consequence
       """.stripMargin

  private val writeDisabled = disabledMessage("Write to Cache", "WILL NOT be cached")
  private val readDisabled = disabledMessage("Read from Cache", "WILL ALL be executed")
}
