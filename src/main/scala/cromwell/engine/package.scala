package cromwell

import java.nio.file.{Path, Paths}
import java.util.UUID

import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue
import cromwell.engine.backend.Backend
import cromwell.engine.workflow.WorkflowOptions
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}
import spray.json._
import lenthall.config.ScalaConfig._

import scala.language.implicitConversions
import scala.util.{Failure, Success, Try}

/**
 * ==Cromwell Execution Engine==
 *
 * Given a WDL file and a backend to execute on, this package provides an API to launch a workflow
 * and query its progress.
 *
 * Internally, this package is built on top of [[cromwell.binding]].
 */
package object engine {

  private val DefaultCallCachingValue = false

  case class WorkflowId(id: UUID) {
    override def toString = id.toString
    def shortString = id.toString.split("-")(0)
  }

  object WorkflowId {
    def fromString(id: String): WorkflowId = new WorkflowId(UUID.fromString(id))
    def randomId() = WorkflowId(UUID.randomUUID())
  }

  object WorkflowDescriptor {
    private def disabledMessage(readWrite: String, consequence: String) =
      s"""$readWrite is enabled in the workflow options but Call Caching is disabled in this Cromwell instance.
         |As a result the calls in this workflow $consequence
       """.stripMargin

    val writeDisabled = disabledMessage("Write to Cache", "WILL NOT be cached")
    val readDisabled = disabledMessage("Read from Cache", "WILL ALL be executed")
  }

  /**
   * Constructs a representation of a particular workflow invocation.  As with other
   * case classes and apply() methods, this will throw an exception if it cannot be
   * created
   */
  case class WorkflowDescriptor(id: WorkflowId, sourceFiles: WorkflowSourceFiles) {
    import WorkflowDescriptor._

    // TODO: Extract this from here (there is no need to reload the configuration for each workflow)
    lazy private [engine] val conf = ConfigFactory.load

    val workflowOptions = Try(sourceFiles.workflowOptionsJson.parseJson) match {
      case Success(options: JsObject) => WorkflowOptions.fromJsonObject(options).get // .get here to purposefully throw the exception
      case Success(other) => throw new Throwable(s"Expecting workflow options to be a JSON object, got $other")
      case Failure(ex) => throw ex
    }

    val props = sys.props
    val workflowLogger = props.get("LOG_MODE") match {
      case Some(x) if x.toUpperCase.contains("SERVER") => makeFileLogger(
        Paths.get(props.getOrElse("LOG_ROOT", ".")),
        s"workflow.$id.log",
        Level.toLevel(props.getOrElse("LOG_LEVEL", "debug"))
      )
      case _ => NOPLogger.NOP_LOGGER
    }

    // Call Caching
    // TODO: Add to lenthall
    def getConfigOption(key: String): Option[Config] = if (conf.hasPath(key)) Option(conf.getConfig(key)) else None
    private lazy val configCallCaching = getConfigOption("call-caching") map { _.getBooleanOr("enabled", DefaultCallCachingValue) } getOrElse DefaultCallCachingValue
    private lazy val optionCacheWriting = workflowOptions.getBoolean("write-to-cache") getOrElse configCallCaching
    private lazy val optionCacheReading = workflowOptions.getBoolean("read-from-cache") getOrElse configCallCaching

    if (!configCallCaching) {
      if (optionCacheWriting) logWriteDisabled()
      if (optionCacheReading) logReadDisabled()
    }

    lazy val writeToCache = configCallCaching && optionCacheWriting
    lazy val readFromCache = configCallCaching && optionCacheReading

    val backend = Backend.from(workflowOptions.getOrElse("default_backend", conf.getConfig("backend").getString("backend")))
    val namespace = NamespaceWithWorkflow.load(sourceFiles.wdlSource, backend.backendType)
    val name = namespace.workflow.unqualifiedName
    val shortId = id.toString.split("-")(0)

    backend.assertWorkflowOptions(workflowOptions)

    val rawInputs = Try(sourceFiles.inputsJson.parseJson) match {
      case Success(JsObject(inputs)) => inputs
      case _ => throw new Throwable(s"Workflow ${id.toString} contains bad inputs JSON: ${sourceFiles.inputsJson}")
    }

    val IOInterface = backend.ioInterface(workflowOptions)

    // Currently we are throwing an exception if construction of the workflow descriptor fails, hence .get on the Trys
    val coercedInputs = namespace.coerceRawInputs(rawInputs).get
    val declarations = namespace.staticDeclarationsRecursive(coercedInputs, backend.engineFunctions(IOInterface)).get
    val actualInputs: WorkflowCoercedInputs = coercedInputs ++ declarations

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

  /**
   * Represents the collection of source files that a user submits to run a workflow
   */
  case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  case class AbortFunction(function: ()=>Unit)
  case class AbortRegistrationFunction(register: AbortFunction=>Unit)

  sealed trait WorkflowState {
    def isTerminal: Boolean
  }

  private lazy val workflowStates = Seq(WorkflowSubmitted, WorkflowRunning, WorkflowFailed, WorkflowSucceeded, WorkflowAborting, WorkflowAborted)

  object WorkflowState {
    def fromString(str: String): WorkflowState = workflowStates.find(_.toString == str).getOrElse(
      throw new NoSuchElementException(s"No such WorkflowState: $str"))
  }

  case object WorkflowSubmitted extends WorkflowState {
    override val toString: String = "Submitted"
    override val isTerminal = false
  }
  
  case object WorkflowRunning extends WorkflowState {
    override val toString: String = "Running"
    override val isTerminal = false
  }

  case object WorkflowAborting extends WorkflowState {
    override val toString: String = "Aborting"
    override val isTerminal = false
  }

  case object WorkflowFailed extends WorkflowState {
    override val toString: String = "Failed"
    override val isTerminal = true
  }
  
  case object WorkflowSucceeded extends WorkflowState {
    override val toString: String = "Succeeded"
    override val isTerminal = true
  }

  case object WorkflowAborted extends WorkflowState {
    override val toString: String = "Aborted"
    override val isTerminal = true
  }

  object SymbolStoreEntry {
    private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
      val lastIndex = fullyQualifiedName.lastIndexOf(".")
      (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
    }

    def apply(fullyQualifiedName: FullyQualifiedName, wdlValue: WdlValue, symbolHash: SymbolHash, input: Boolean): SymbolStoreEntry = {
      val (scope, name) = splitFqn(fullyQualifiedName)
      val key = SymbolStoreKey(scope, name, index = None, input)
      SymbolStoreEntry(key, wdlValue.wdlType, Option(wdlValue), Option(symbolHash))
    }

    def toWorkflowOutputs(t: Traversable[SymbolStoreEntry]): WorkflowOutputs = t.map { e =>
      s"${e.key.scope}.${e.key.name}" -> CallOutput(e.wdlValue.get, e.symbolHash.get)
    }.toMap

    def toCallOutputs(traversable: Traversable[SymbolStoreEntry]): CallOutputs = traversable.map { entry =>
      entry.key.name -> CallOutput(entry.wdlValue.get, entry.symbolHash.get)
    }.toMap
  }

  case class SymbolStoreKey(scope: String, name: String, index: Option[Int], input: Boolean) {
    def fqn: FullyQualifiedName = s"$scope.$name"
  }

  case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue], symbolHash: Option[SymbolHash]) {
    def isInput: Boolean = key.input
    def isOutput: Boolean = !isInput
    def scope: String = key.scope
  }

  object ExecutionStatus extends Enumeration {
    type ExecutionStatus = Value
    val NotStarted, Starting, Running, Failed, Done, Aborted = Value

    implicit class EnhancedExecutionStatus(val status: ExecutionStatus) extends AnyVal {
      def isTerminal: Boolean = {
        Seq(Failed, Done, Aborted) contains status
      }
    }

    implicit class EnhancedString(val string: String) extends AnyVal {
      def toExecutionStatus: ExecutionStatus = ExecutionStatus.withName(string)
    }
  }

  /*
   * Type and implicit conversion classes for ExecutionIndex
   */
  object ExecutionIndex {
    type ExecutionIndex = Option[Int]
    val IndexNone = -1 // "It's a feature" https://bugs.mysql.com/bug.php?id=8173

    implicit class IndexEnhancedInt(val value: Int) extends AnyVal {
      def toIndex: ExecutionIndex = value match {
        case IndexNone => None
        case i => Option(i)
      }
    }

    implicit class IndexEnhancedIndex(val index: ExecutionIndex) extends AnyVal {
      def fromIndex: Int = index match {
        case None => IndexNone
        case Some(i) => i
      }
    }

    implicit val ExecutionIndexOrdering = new Ordering[ExecutionIndex] {
      override def compare(x: ExecutionIndex, y: ExecutionIndex): Int = {
        x.fromIndex.compareTo(y.fromIndex)
      }
    }
  }
}
