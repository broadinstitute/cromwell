package cromwell.engine

import java.nio.file._

import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.Backend
import cromwell.engine.backend.DefaultWorkflowEngineFunctions
import cromwell.engine.db.DataAccess.globalDataAccess
import cromwell.engine.io.gcs.GoogleCloudStorage
import cromwell.engine.io.shared.SharedFileSystemIoInterface
import cromwell.engine.io.{IoInterface, IoManager}
import cromwell.engine.workflow.WorkflowOptions
import cromwell.logging.WorkflowLogger
import cromwell.util.{SimpleExponentialBackoff, TryUtil}
import lenthall.config.ScalaConfig._
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}
import spray.json.{JsObject, _}
import wdl4s._
import wdl4s.values.{WdlFile, WdlSingleFile}

import scala.concurrent._
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._

case class WorkflowDescriptor(id: WorkflowId,
                              sourceFiles: Option[WorkflowSourceFiles],
                              workflowOptions: WorkflowOptions,
                              namespace: NamespaceWithWorkflow,
                              coercedInputs: WorkflowCoercedInputs,
                              declarations: WorkflowCoercedInputs,
                              configCallCaching: Boolean,
                              lookupDockerHash: Boolean,
                              ioManager: IoInterface,
                              wfContext: WorkflowContext,
                              engineFunctions: WorkflowEngineFunctions) {
  import WorkflowDescriptor._

  val shortId = id.toString.split("-")(0)
  val name = namespace.workflow.unqualifiedName
  val actualInputs: WorkflowCoercedInputs = coercedInputs ++ declarations
  val props = sys.props
  private val relativeWorkflowRootPath = s"$name/$id"
  private val log = WorkflowLogger("WorkflowDescriptor", this)
  val workflowOutputsPath = workflowOptions.get("outputs_path") recover { case e: IllegalArgumentException =>
    log.warn("outputs_path expected to be of type String", e)
    throw e
  }
  lazy val fileHasher: FileHasher = { wdlFile: WdlFile => SymbolHash(ioManager.hash(wdlFile.value)) }
  private lazy val optionCacheWriting = workflowOptions.getBoolean("write_to_cache") getOrElse configCallCaching
  private lazy val optionCacheReading = workflowOptions.getBoolean("read_from_cache") getOrElse configCallCaching

  if (!configCallCaching) {
    if (optionCacheWriting) logWriteDisabled()
    if (optionCacheReading) logReadDisabled()
  }

  lazy val writeToCache = configCallCaching && optionCacheWriting
  lazy val readFromCache = configCallCaching && optionCacheReading

  lazy val workflowLogger = props.get("LOG_MODE") match {
    case Some(x) if x.toUpperCase.contains("SERVER") => makeFileLogger(
      Paths.get(props.getOrElse("LOG_ROOT", ".")),
      s"workflow.$id.log",
      Level.toLevel(props.getOrElse("LOG_LEVEL", "debug"))
    )
    case _ => NOPLogger.NOP_LOGGER
  }

  lazy val workflowRootPath = wfContext.path
  def workflowRootPathWithBaseRoot(rootPath: String): Path = Paths.get(WorkflowDescriptor.buildWorkflowRootPath(rootPath, name, id))

  def copyWorkflowOutputs(implicit executionContext: ExecutionContext): Future[Unit] = {
    // Try to copy outputs to final destination
    workflowOutputsPath map copyOutputFiles getOrElse Future.successful(())
  }

  private def copyOutputFiles(destDirectory: String)(implicit executionContext: ExecutionContext): Future[Unit] = {

    def copyFile(file: WdlFile): Try[Unit] = {
      val src = Paths.get(file.valueString)
      val wfPath = wfContext.path.toAbsolutePath
      val relativeFilePath = Paths.get(relativeWorkflowRootPath).resolve(src.subpath(wfPath.getNameCount, src.getNameCount))
      val dest = Paths.get(destDirectory).resolve(relativeFilePath)

      def copy(): Unit = {
        log.info(s"Trying to copy output file $src to $dest")
        Files.createDirectories(dest.getParent)
        Files.copy(src, dest)
      }

      TryUtil.retryBlock(
        fn = (retries: Option[Unit]) => copy(),
        retryLimit = Option(5),
        backoff = SimpleExponentialBackoff(5 seconds, 10 seconds, 1.1D),
        logger = log,
        failMessage = Option(s"Failed to copy file $src to $dest"),
        isFatal = (t: Throwable) => t.isInstanceOf[FileAlreadyExistsException]
      ) recover {
        case _: FileAlreadyExistsException => log.info(s"Tried to copy the same file multiple times. Skipping subsequent copies for $src")
      }
    }

    def processOutputs(outputs: Traversable[SymbolStoreEntry]): Unit = {
      // All outputs should have wdl values at this point, if they don't there's nothing we can do here
      val copies = TryUtil.sequence(outputs map { o => Try(o.wdlValue.get) } toSeq) match {
        case Success(wdlValues) => wdlValues flatMap { _ collectAsSeq { case f: WdlSingleFile => f } } map copyFile
        case Failure(e) => throw new Throwable(s"Unable to resolve the following workflow outputs for workflow $id: ${e.getMessage}")
      }

      // Throw an exception if one or more of the copies failed.
      TryUtil.sequence(copies) match {
        case Success(_) => ()
        case Failure(e) => throw new Throwable(s"Output copy failed for the following files:\n ${e.getMessage}")
      }
    }

    globalDataAccess.getWorkflowOutputs(id) map processOutputs
  }

  private def makeFileLogger(root: Path, name: String, level: Level): Logger = {
    val ctx = LoggerFactory.getILoggerFactory.asInstanceOf[LoggerContext]
    val encoder = new PatternLayoutEncoder()
    encoder.setPattern("%date %-5level - %msg%n")
    encoder.setContext(ctx)
    encoder.start()

    val path = root.resolve(name).toAbsolutePath.toString
    val appender = new FileAppender[ILoggingEvent]()
    appender.setFile(path)
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

  private def logWriteDisabled() = workflowLogger.warn(writeDisabled)
  private def logReadDisabled() = workflowLogger.warn(readDisabled)
}

object WorkflowDescriptor {
  //TODO: This constructor is there because of dependencies in SlickDataAccess. Should use just one construction (i.e. the second one), where it is assumed the WF passed validations
  def apply(id: WorkflowId, sourceFiles: WorkflowSourceFiles, conf: Config): WorkflowDescriptor = {
    val namespace = NamespaceWithWorkflow.load(sourceFiles.wdlSource)
    val inputs = sourceFiles.inputsJson.parseJson.asJsObject.fields
    val coercedInputs = namespace.coerceRawInputs(inputs).getOrElse(throw new IllegalArgumentException(s"Failed to coerce inoputs for WF $id"))
    WorkflowDescriptor(id, Option(sourceFiles), namespace, Option(coercedInputs), Option(sourceFiles.workflowOptionsJson), conf)
  }

  def apply(id: WorkflowId,
            sourceFiles: Option[WorkflowSourceFiles],
            namespaceWithWorkflow: NamespaceWithWorkflow,
            coercedInputs: Option[WorkflowCoercedInputs],
            workflowOptionsJson: Option[WorkflowOptionsJson],
            conf: Config = ConfigFactory.load()): WorkflowDescriptor = {
    val wfContext = new WorkflowContext(id.toString)
    val ioManager = new IoManager(Seq(Option(SharedFileSystemIoInterface.instance)).flatten)
    val engineFunctions = new DefaultWorkflowEngineFunctions(ioManager, wfContext)
    val wfOptions = validateWorkflowOptions(id, workflowOptionsJson.getOrElse("{}"))
    val decl = validateDeclarations(id, namespaceWithWorkflow, coercedInputs.get, engineFunctions)
    new WorkflowDescriptor(id, sourceFiles, wfOptions, namespaceWithWorkflow, coercedInputs.get, decl, configCallCaching(conf), lookupDockerHash(conf),ioManager, wfContext, engineFunctions)
  }

  private def buildWorkflowRootPath(rootPath: String, name: String, workflowId: WorkflowId) = s"$rootPath/$name/$workflowId"

  //TODO: Move it to the ValidateActor part, can't move now since this has engine functions (i.e. Backend) dependency.
  private def validateDeclarations(id: WorkflowId,
                                   namespace: NamespaceWithWorkflow,
                                   coercedInputs: WorkflowCoercedInputs,
                                   engineFunctions: WorkflowEngineFunctions): WorkflowCoercedInputs = {
    namespace.staticDeclarationsRecursive(coercedInputs, engineFunctions) match {
      case Success(d) => d
      case Failure(e) => throw new IllegalStateException(s"Workflow $id has invalid declarations: ${e.getMessage}")
    }
  }

  //TODO: Same as above, should be a part of individual backends
  private def validateWorkflowOptions(id: WorkflowId,
                                      optionsJson: WorkflowOptionsJson): WorkflowOptions = {
    WorkflowOptions.fromJsonString(optionsJson) match {
      case Success(o) => o
      case Failure(e) =>
        throw new IllegalArgumentException(s"Workflow $id contains bad options JSON: ${e.getMessage}", e)
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
}
