package cromwell

import java.util.UUID

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue
import cromwell.engine.backend.{ExecutionHandle, ExecutionResult}
import org.joda.time.DateTime

import scala.concurrent.{ExecutionContext, Future}
import scala.language.implicitConversions
import scalaz.ValidationNel

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
  private val DefaultLookupDockerHash = false

  case class WorkflowId(id: UUID) {
    override def toString = id.toString
    def shortString = id.toString.split("-")(0)
  }

  object WorkflowId {
    def fromString(id: String): WorkflowId = new WorkflowId(UUID.fromString(id))
    def randomId() = WorkflowId(UUID.randomUUID())
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

  case class ExecutionEventEntry(description: String, startTime: DateTime, endTime: DateTime)

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

  type ErrorOr[+A] = ValidationNel[String, A]

  implicit class EnhancedFutureFuture[A](val ffa: Future[Future[A]])(implicit ec: ExecutionContext) {
    def flatten: Future[A] = ffa flatMap { fa => fa }
  }

  implicit class EnhancedExecutionHandle(val handle: ExecutionHandle) extends AnyVal {
    def future = Future.successful(handle)
  }

  implicit class EnhancedExecutionResult(val result: ExecutionResult) extends AnyVal {
    def future = Future.successful(result)
  }

  final case class ExecutionHash(overallHash: String, dockerHash: Option[String])
}
