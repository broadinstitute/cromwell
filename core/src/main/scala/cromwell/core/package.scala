package cromwell

import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s._
import wdl4s.expression.WdlFunctions
import wdl4s.types.WdlType
import wdl4s.values.{WdlValue, SymbolHash}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz._

package object core {
  // root can be a Path instead of a String in PBE. stdout / err too but it doesn't really bring values since they're just stringified to WdlFiles
  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  // This could move to WDL4s
  type StringMapper = String => String
  type WdlValueMapper = WdlValue => WdlValue

  case class ResolvedDeclaration(wdlType: WdlType, name: String, unEvaluatedValue: WdlValue) {
    def toNameValuePair = (name, unEvaluatedValue)
  }

  // Can evaluate a wdl value
  class Evaluator(defaultLookup: ScopedLookupFunction,
                  engineFunctions: WdlFunctions[WdlValue],
                  preValueMapper: StringMapper,
                  postValueMapper: WdlValueMapper) {

    final private case class AttemptedLookupResult(name: String, value: Try[WdlValue]) { def toPair = name -> value }

    private def lookup = defaultLookup compose preValueMapper

    // Evaluate a wdl value with the given lookup function, and optionally coerce the result to coerceTo
    private def evaluateValueWith(customLookup: ScopedLookupFunction, wdlValue: WdlValue) = {
      wdlValue match {
        case wdlExpression: WdlExpression => wdlExpression.evaluate(customLookup, engineFunctions) map postValueMapper
        case v: WdlValue => Success(v)
      }
    }

    private def evaluateDeclarationWith(customLookup: ScopedLookupFunction, declaration: ResolvedDeclaration) = {
      evaluateValueWith(customLookup, declaration.unEvaluatedValue) flatMap declaration.wdlType.coerceRawValue
    }

    private def mapLookup(values: Map[LocallyQualifiedName, Try[WdlValue]], identifier: String): Try[WdlValue] = {
      values.getOrElse(identifier, Failure(new WdlExpressionException(s"Could not resolve variable '$identifier' as an input parameter")))
    }

    def evaluateValue(wdlValue: WdlValue): Try[WdlValue] = evaluateValueWith(lookup, wdlValue)

    def evaluateDeclaration(declaration: ResolvedDeclaration) = evaluateDeclarationWith(lookup, declaration)

    def evaluateDeclarations(declarations: Seq[ResolvedDeclaration]): Map[LocallyQualifiedName, Try[WdlValue]] = {

      declarations.foldLeft(Map.empty[LocallyQualifiedName, Try[WdlValue]])((evaluatedValues, declaration) => {

        def enhancedLookup(identifier: String) = mapLookup(evaluatedValues, identifier) getOrElse lookup(identifier)

        evaluatedValues + (declaration.name -> evaluateDeclarationWith(enhancedLookup, declaration))
      })
    }

    def evaluateValues(wdlValues: Map[LocallyQualifiedName, WdlValue]): Map[LocallyQualifiedName, Try[WdlValue]] = {

      wdlValues.foldLeft(Map.empty[LocallyQualifiedName, Try[WdlValue]])((evaluatedValues, pair) => {
        val (name, wdlValue) = pair

        def enhancedLookup(identifier: String) = mapLookup(evaluatedValues, identifier) getOrElse lookup(identifier)

       evaluatedValues + (name -> evaluateValueWith(enhancedLookup, wdlValue))
      })
    }
  }

  // Can build an evaluator from engine functions and valueMapper
  class EvaluatorBuilder(builder: (WdlFunctions[WdlValue], StringMapper, WdlValueMapper) => Evaluator) {

    def build(functions: WdlFunctions[WdlValue],
              preMapper: StringMapper = identity[String],
              postMapper: WdlValueMapper = identity[WdlValue]) = builder(functions, preMapper, postMapper)
  }

  type ErrorOr[+A] = ValidationNel[String, A]
  type LocallyQualifiedName = String
  case class CallOutput(wdlValue: WdlValue, hash: Option[SymbolHash])
  type CallOutputs = Map[LocallyQualifiedName, CallOutput]
  type EvaluatedRuntimeAttributes = Map[String, WdlValue]

  implicit class EnhancedFQN(val fqn: FullyQualifiedName) extends AnyVal {
    def unqualified: LocallyQualifiedName = splitFqn._2

    def splitFqn: (String, String) = {
      val lastIndex = fqn.lastIndexOf(".")
      if (lastIndex == -1) ("", fqn)
      else (fqn.substring(0, lastIndex), fqn.substring(lastIndex + 1))
    }
  }
}