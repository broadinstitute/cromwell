package wom

import org.apache.commons.codec.digest.DigestUtils
import wom.callable.Callable.InputDefinition
import wom.values.{WomFile, WomValue}

import scala.util.Try

trait TsvSerializable {
  def tsvSerialize: Try[String]
}

class WomExpressionException(message: String = null, cause: Throwable = null) extends RuntimeException(message, cause)
final case class OptionalNotSuppliedException(operationName: String) extends Exception(s"Sorry! Operation $operationName is not supported on empty optional values. You might resolve this using select_first([optional, default]) to guarantee that you have a filled value.")

case class JobOutput(womValue: WomValue)

package object values {

  type FileHasher = WomFile => SymbolHash

  type WomEvaluatedCallInputs = Map[InputDefinition, WomValue]

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
  type HostInputs = Map[String, WomValue]
  type EvaluatedRuntimeAttributes = Map[String, WomValue]
  // This one really does seem WOM specific
  type ExecutableInputMap = Map[String, Any]
}

/**
  * @param commandString The string representing the instantiation of this command.
  * @param createdFiles Any files created as side effects of instantiating the command.
  */
case class InstantiatedCommand(commandString: String, createdFiles: List[WomFile] = List.empty)
