package cromwell.engine

import java.nio.file.{Path, Paths}

import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import com.typesafe.config.{Config, ConfigFactory}
import wdl4s._
import wdl4s.values.WdlFile
import cromwell.engine.backend.{DefaultWorkflowEngineFunctions, BackendType, Backend}
import cromwell.engine.io.{IoInterface, IoManager}
import cromwell.engine.io.gcs.GoogleCloudStorage
import cromwell.engine.io.shared.SharedFileSystemIoInterface
import cromwell.engine.workflow.WorkflowOptions
import lenthall.config.ScalaConfig._
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}
import spray.json.{JsObject, _}

import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._

case class WorkflowDescriptor(id: WorkflowId,
                              sourceFiles: WorkflowSourceFiles,
                              workflowOptions: WorkflowOptions,
                              rawInputs: Map[String, JsValue],
                              namespace: NamespaceWithWorkflow,
                              coercedInputs: WorkflowCoercedInputs,
                              declarations: WorkflowCoercedInputs,
                              configCallCaching: Boolean,
                              lookupDockerHash: Boolean,
                              gcsInterface: Try[GoogleCloudStorage],
                              ioManager: IoInterface,
                              wfContext: WorkflowContext,
                              engineFunctions: WorkflowEngineFunctions) {
  import WorkflowDescriptor._

  val shortId = id.toString.split("-")(0)
  val name = namespace.workflow.unqualifiedName
  val actualInputs: WorkflowCoercedInputs = coercedInputs ++ declarations
  val props = sys.props
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
  def apply(id: WorkflowId, sourceFiles: WorkflowSourceFiles): WorkflowDescriptor = {
    WorkflowDescriptor(id, sourceFiles, ConfigFactory.load)
  }

  def apply(id: WorkflowId, sourceFiles: WorkflowSourceFiles, conf: Config): WorkflowDescriptor = {
    validateWorkflowDescriptor(id, sourceFiles, conf) match {
      case scalaz.Success(w) => w
      case scalaz.Failure(f) =>
        throw new IllegalArgumentException(s"""Workflow $id failed to process inputs:\n${f.toList.mkString("\n")}""")
    }
  }

  private def validateWorkflowDescriptor(id: WorkflowId,
                                         sourceFiles: WorkflowSourceFiles,
                                         conf: Config): ErrorOr[WorkflowDescriptor] = {
    val namespace = validateNamespace(id, sourceFiles.wdlSource)
    val rawInputs = validateRawInputs(id, sourceFiles.inputsJson)
    val options = validateWorkflowOptions(id, sourceFiles.workflowOptionsJson)

    val runtimeAttributes = for {
      n <- namespace.disjunction
    } yield validateRuntimeAttributes(id, n)

    (namespace |@| rawInputs |@| options |@| runtimeAttributes.validation) { (_, _, _, _) } match {
      case scalaz.Success((n, r, o, a)) => buildWorkflowDescriptor(id, sourceFiles, n, r, o, conf)
      case scalaz.Failure(f) => f.toList.mkString("\n").failureNel
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFiles,
                                      namespace: NamespaceWithWorkflow,
                                      rawInputs: Map[String, JsValue],
                                      options: WorkflowOptions,
                                      conf: Config): ErrorOr[WorkflowDescriptor] = {
    val gcsInterface = GoogleCloudStorage.userAuthenticated(options) orElse GoogleCloudStorage.cromwellAuthenticated
    val ioManager = new IoManager(Seq(gcsInterface.toOption, Option(SharedFileSystemIoInterface.instance)).flatten)
    val wfContext = new WorkflowContext(id.toString)
    val engineFunctions = new DefaultWorkflowEngineFunctions(ioManager, wfContext)

    val validatedDescriptor = for {
      c <- validateCoercedInputs(id, rawInputs, namespace).disjunction
      d <- validateDeclarations(id, namespace, options, c, engineFunctions).disjunction
    } yield WorkflowDescriptor(id, sourceFiles, options, rawInputs, namespace, c, d, configCallCaching(conf), lookupDockerHash(conf),
      gcsInterface, ioManager, wfContext, engineFunctions)
    validatedDescriptor.validation
  }

  private def validateNamespace(id: WorkflowId, source: WdlSource): ErrorOr[NamespaceWithWorkflow] = {
    try {
      NamespaceWithWorkflow.load(source).successNel
    } catch {
      case e: Exception => s"Workflow $id unable to load namespace: ${e.getMessage}".failureNel
    }
  }

  private def validateRuntimeAttributes(id: WorkflowId, namespace: NamespaceWithWorkflow): ErrorOr[Unit] = {
    Try(namespace.workflow.calls foreach { x => x.task.runtimeAttributes }) match {
      case scala.util.Success(_) => ().successNel
      case scala.util.Failure(e) => s"Workflow $id contains bad runtime attributes: ${e.getMessage}".failureNel
    }
  }

  private def validateWorkflowOptions(id: WorkflowId,
                                      optionsJson: WorkflowOptionsJson): ErrorOr[WorkflowOptions] = {
    WorkflowOptions.fromJsonString(optionsJson) match {
      case Success(o) => validateBackendOptions(id, o)
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
      case Failure(e) => e.getMessage.failureNel
    }
  }

  private def validateBackendOptions(id: WorkflowId, options: WorkflowOptions): ErrorOr[WorkflowOptions] = {
    try {
      options.successNel
    } catch {
      case e: Exception => s"Workflow $id has invalid options: ${e.getMessage}".failureNel
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
}