package wom

import org.apache.commons.codec.digest.DigestUtils
import wom.callable.Callable.InputDefinition
import wom.graph.LocalName
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
  
  implicit class EnhancedWomEvaluatedCallInputs(val inputs: WomEvaluatedCallInputs) extends AnyVal {
    def prettyString = inputs.map({
      case (inputDef, womValue) => s"${inputDef.name} -> ${womValue.valueString}"
    }).mkString(", ")
  }

  implicit class HashableString(val value: String) extends AnyVal with Hashable {
    def md5Sum: String = DigestUtils.md5Hex(value)
    def md5SumShort: String = value.md5Sum.substring(0, 8)
  }

  implicit class ShellQuoteHelper(val s: String) extends AnyVal {
    // Escape any double quote characters in the string and surround with double quotes. The CWL spec doesn't
    // appear to say whether quotes should be of the single or double variety, but the conformance tests clearly
    // expect double quotes to allow for interpolation of environment variables.
    def shellQuote: String = '"' + s.replaceAll("\"", "\\\"") + '"'
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
  * @param environmentVariables Key/value environment variable pairs.
  * @param createdFiles Any files created as side effects of instantiating the command.
  * @param stdinRedirection An optional redirection of standard input from a stringified filename.
  * @param preprocessedInputs A List of tuples of preprocessed inputs. This is a List rather than a Map because this class
  *                   needs to have a monoid instance. If this class contained Maps there would need to be monoid
  *                   instances for their WomValue values which is unpossible (and the way Cromwell uses this data
  *                   is also not actually necessary).
  * @param valueMappedPreprocessedInputs A List of tuples of preprocessed and value-mapped task inputs.
  */
final case class InstantiatedCommand(commandString: String,
                                     environmentVariables: Map[String, String] = Map.empty,
                                     createdFiles: List[CommandSetupSideEffectFile] = List.empty,
                                     stdinRedirection: Option[String] = None,
                                     preprocessedInputs: List[(LocalName, WomValue)] = List.empty,
                                     valueMappedPreprocessedInputs: List[(LocalName, WomValue)] = List.empty)

/**
  * File created as a side effect of instantiating the command.
  *
  * @param file The WomFile which we need to localize
  * @param relativeLocalPath Optionally, an override of the usual localization logic. A path relative to execution root.
  */
final case class CommandSetupSideEffectFile(file: WomFile, relativeLocalPath: Option[String] = None)
