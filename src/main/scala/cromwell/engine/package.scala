package cromwell

import wdl4s._
import org.joda.time.DateTime
import wdl4s.values.{WdlValue, WdlFile}

import scala.language.implicitConversions
import scalaz.ValidationNel

package object engine {

  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  /**
   * Represents the collection of source files that a user submits to run a workflow
   */
  final case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  final case class AbortFunction(function: ()=>Unit)
  final case class AbortRegistrationFunction(register: AbortFunction=>Unit)

  final case class ExecutionEventEntry(description: String, startTime: DateTime, endTime: DateTime)

  type ErrorOr[+A] = ValidationNel[String, A]

  case class SymbolHash(value: String) extends Ordered[SymbolHash] {
    def compare(that: SymbolHash) = this.value compare that.value
  }
  type FileHasher = WdlFile => SymbolHash

  type WorkflowOptionsJson = String
  type WorkflowOutputs = Map[FullyQualifiedName, CallOutput]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type CallInputs = Map[String, WdlValue]
  case class CallOutput(wdlValue: WdlValue, hash: Option[SymbolHash])
  type CallOutputs = Map[LocallyQualifiedName, CallOutput]
  type HostInputs = Map[String, WdlValue]

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def isScatter = fqn.contains(Scatter.FQNIdentifier)
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, CallOutput]) extends AnyVal {
    def mapToValues: Map[A, WdlValue] = m map {
      case (k, CallOutput(wdlValue, hash)) => (k, wdlValue)
    }
  }
}
