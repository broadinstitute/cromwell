package cromwell

import cromwell.binding.values.{WdlFile, WdlValue}

import scala.language.implicitConversions

/**
 * ==WDL Bindings for Scala==
 *
 * This package contains an API to convert WDL ASTs into native Scala objects.
 *
 * Clients of this package should be able to interface with a WDL file without ever coming in direct contact with an Abstract Syntax Tree
 *
 * Internally, this package is built on top of [[cromwell.parser]].
 */

package object binding {
  type WdlSource = String
  type WdlJson = String
  type WorkflowOptionsJson = String
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type WorkflowOutputs = Map[FullyQualifiedName, CallOutput]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type CallInputs = Map[String, WdlValue]
  case class CallOutput(wdlValue: WdlValue, hash: SymbolHash)
  type CallOutputs = Map[LocallyQualifiedName, CallOutput]
  type HostInputs = Map[String, WdlValue]
  type ImportResolver = String => WdlSource
  case class SymbolHash(value: String) extends Ordered[SymbolHash] {
    def compare(that: SymbolHash) = this.value compare that.value
  }
  type FileHasher = WdlFile => SymbolHash

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
  
  object Patterns {

    val WorkflowName = """
        (?x)                             # Turn on comments and whitespace insensitivity.

        (                                # Begin capture.

          [a-zA-Z][a-zA-Z0-9_]*          # WDL identifier naming pattern of an initial alpha character followed by zero
                                         # or more alphanumeric or underscore characters.

        )                                # End capture.
      """.trim.r

    val CallFullyQualifiedName = """
      (?x)                               # Turn on comments and whitespace insensitivity.

      (                                  # Begin outer capturing group for FQN.

        (?:[a-zA-Z][a-zA-Z0-9_]*)        #   Inner noncapturing group for top-level workflow name. This is the WDL
                                         #   identifier naming pattern of an initial alpha character followed by zero
                                         #   or more alphanumeric or underscore characters.

        (?:\.[a-zA-Z][a-zA-Z0-9_]*){1}   #   Inner noncapturing group for call name, a literal dot followed by a WDL
                                         #   identifier.  Currently this is quantified to {1} since the call name is
                                         #   mandatory and nested workflows are not supported.  This could be changed
                                         #   to + or a different quantifier if these assumptions change.

      )                                  # End outer capturing group for FQN.


      (?:                                # Begin outer noncapturing group for shard.
        \.                               #   Literal dot.
        (\d+)                            #   Captured shard digits.
      )?                                 # End outer optional noncapturing group for shard.
      """.trim.r                         // The trim is necessary as (?x) must be at the beginning of the regex.
  }
}
