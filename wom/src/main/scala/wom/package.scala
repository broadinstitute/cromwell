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
    // `cwltool` invokes `shellescape.quote` to escape its arguments. This should reimplement that algorithm.
    // [1] https://github.com/common-workflow-language/cwltool/blob/93a122a3da7926be8efe9a0952214ee7244d0d46/cwltool/draft2tool.py#L495-L496
    // [2] https://github.com/chrissimpkins/shellescape/blob/1600bceca32a8df2d8928751003f429dad9fe77e/lib/shellescape/main.py#L21
    def shellQuote: String = "'" + s.replaceAll("'", "'\"'\"'") + "'"
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
  type WorkflowUrl = String
  type WorkflowTypeVersion = String
  type HostInputs = Map[String, WomValue]
  type EvaluatedRuntimeAttributes = Map[String, WomValue]
  // This one really does seem WOM specific
  type ExecutableInputMap = Map[String, Any]
}

/***
  * @param importPath The string representing the resolved import.
  */
final case class ResolvedImportRecord(importPath: String) extends AnyVal

/**
  * @param commandString The string representing the instantiation of this command.
  * @param environmentVariables Key/value environment variable pairs.
  * @param createdFiles Any files created as side effects of instantiating the command.
  * @param evaluatedStdinRedirection An optional redirection of standard input from a stringified filename.
  * @param evaluatedStdoutOverride An optional override of standard output to a specified filename. Standard output will
  *                       always be redirected to a file, this parameter only controls the name of the file.
  * @param evaluatedStderrOverride An optional override of standard error to a specified filename. Standard error will
  *                       always be redirected to a file, this parameter only controls the name of the file.
  * @param preprocessedInputs A List of tuples of preprocessed inputs. This is a List rather than a Map because this class
  *                   needs to have a monoid instance. If this class contained Maps there would need to be monoid
  *                   instances for their WomValue values which is unpossible (and the way Cromwell uses this data
  *                   is also not actually necessary).
  * @param valueMappedPreprocessedInputs A List of tuples of preprocessed and value-mapped task inputs.
  */
final case class InstantiatedCommand(commandString: String,
                                     environmentVariables: Map[String, String] = Map.empty,
                                     createdFiles: List[CommandSetupSideEffectFile] = List.empty,
                                     evaluatedStdinRedirection: Option[String] = None,
                                     evaluatedStdoutOverride: Option[String] = None,
                                     evaluatedStderrOverride: Option[String] = None,
                                     preprocessedInputs: List[(LocalName, WomValue)] = List.empty,
                                     valueMappedPreprocessedInputs: List[(LocalName, WomValue)] = List.empty)

/**
  * File created as a side effect of instantiating the command.
  *
  * @param file The WomFile which we need to localize
  * @param relativeLocalPath Optionally, an override of the usual localization logic. A path relative to execution root.
  */
final case class CommandSetupSideEffectFile(file: WomFile, relativeLocalPath: Option[String] = None)

// Reference to program source code. Allows:
// 1) Relating a WOM error back to the original source code.
// 2) Getting back information lost during compilation.
// For example, in a program like:
//  workflow {
//    call A
//    call B
//    call C
//  }
// the WOM graph does not say which call (A,B, or C) comes first.
//
//
// We would actually like to have the entire source extent covered
// by the AST. However, it is tricky to get full information from Hermes
// at the moment.
final case class SourceFileLocation(startLine: Int)
