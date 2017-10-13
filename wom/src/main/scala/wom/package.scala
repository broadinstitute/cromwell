package wom

import org.apache.commons.codec.digest.DigestUtils
import wom.callable.Callable.InputDefinition
import wom.values.WdlValue

import scala.util.Try

trait TsvSerializable {
  def tsvSerialize: Try[String]
}

class WdlExpressionException(message: String = null, cause: Throwable = null) extends RuntimeException(message, cause)
final case class OptionalNotSuppliedException(operationName: String) extends Exception(s"Sorry! Operation $operationName is not supported on empty optional values. You might resolve this using select_first([optional, default]) to guarantee that you have a filled value.")

case class JobOutput(wdlValue: WdlValue)

package object values {

  type FileHasher = WdlFile => SymbolHash

  type WomEvaluatedCallInputs = Map[InputDefinition, WdlValue]

  implicit class HashableString(val value: String) extends AnyVal with Hashable {
    def md5Sum: String = DigestUtils.md5Hex(value)
    def md5SumShort: String = value.md5Sum.substring(0, 8)
  }
}

package object core {
  type LocallyQualifiedName = String
  type FullyQualifiedName = String
  type WorkflowJson = String
  type WorkflowOutputs = Map[FullyQualifiedName, JobOutput]
  type WorkflowOptionsJson = String
  type WorkflowType = String
  type WorkflowSource = String
  type WorkflowTypeVersion = String
  type CallOutputs = Map[String, JobOutput]
  type HostInputs = Map[String, WdlValue]
  type EvaluatedRuntimeAttributes = Map[String, WdlValue]
  // This one really does seem WDL specific
  type ExecutableInputMap = Map[String, Any]
}
