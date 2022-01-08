package cwl

import cats.syntax.either._
import cats.syntax.validated._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values._
import wom.{CommandPart, InstantiatedCommand}

import scala.language.postfixOps

case class CwlExpressionCommandPart(expr: Expression)(hasShellCommandRequirement: Boolean, expressionLib: ExpressionLib) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
    val stringKeyMap = inputsMap map {
      case (LocalName(localName), value) => localName -> valueMapper(value)
    }
    val parameterContext = ParameterContext(
      functions,
      expressionLib,
      inputs = stringKeyMap, runtimeOption = Option(runtimeEnvironment)
    )
    ExpressionEvaluator.eval(expr, parameterContext) map { womValue =>
      List(InstantiatedCommand(valueMapper(womValue).valueString.shellQuote))
    }
  }
}

/**
  * Generates command parts from a CWL Binding
  */
abstract class CommandLineBindingCommandPart(commandLineBinding: CommandLineBinding)(hasShellCommandRequirement: Boolean, expressionLib: ExpressionLib) extends CommandPart {


  private lazy val prefixAsString = commandLineBinding.prefix.getOrElse("")
  private lazy val prefixAsList = commandLineBinding.prefix.toList
  private lazy val separate = commandLineBinding.effectiveSeparate

  def handlePrefix(value: String) = {
    if (separate) prefixAsList :+ value else List(s"$prefixAsString$value")
  }

  /**
    * Value bound to this command part.
    * InputCommandLineBindingCommandPart has one
    * ArgumentCommandLineBindingCommandPart does not
    */
  def boundValue: Option[WomValue]

  // If the bound value is defined but contains an empty optional value, we should not evaluate the valueFrom
  // Conformance test "stage-unprovided-file" tests this behavior
  private lazy val evaluateValueFrom = boundValue.forall {
    case WomOptionalValue(_, None) => false
    case _ => true
  }

  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
      val stringInputsMap = inputsMap map {
        case (LocalName(localName), value) => localName -> valueMapper(value)
      }
    val parameterContext = ParameterContext(
      functions,
      expressionLib,
      inputs = stringInputsMap,
      runtimeOption = Option(runtimeEnvironment),
      self = boundValue.getOrElse(ParameterContext.EmptySelf)
    )

    val evaluatedValueFrom = commandLineBinding.optionalValueFrom flatMap {
      case StringOrExpression.Expression(expression) if evaluateValueFrom => Option(ExpressionEvaluator.eval(expression, parameterContext) map valueMapper)
      case StringOrExpression.String(string) if evaluateValueFrom => Option(WomString(string).validNel)
      case _ => None
    }

    val evaluatedWomValue: Checked[WomValue] = evaluatedValueFrom.orElse(boundValue.map(_.validNel)) match {
      case Some(womValue) => womValue.map(valueMapper).toEither
      case None => "Command line binding has no valueFrom field and no bound value".invalidNelCheck
    }

    def applyShellQuote(value: String): String = commandLineBinding.shellQuote match {
      // Only honor shellQuote = false if ShellCommandRequirement is enabled.
      // Conformance test "Test that shell directives are not interpreted."
      case Some(false) if hasShellCommandRequirement => value
      case _ => value.shellQuote
    }

    def processValue(womValue: WomValue): List[String] = womValue match {
      case WomOptionalValue(_, Some(value)) => processValue(valueMapper(value))
      case WomOptionalValue(_, None) => List.empty
      case _: WomString | _: WomInteger | _: WomFile | _: WomLong | _: WomFloat => handlePrefix(valueMapper(womValue).valueString)
      // For boolean values, use the value of the boolean to choose whether to print the prefix or not
      case WomBoolean(false) => List.empty
      case WomBoolean(true) => prefixAsList
      case WomArray(_, values) => commandLineBinding.itemSeparator match {
        case Some(_) if values.isEmpty => List.empty
        case Some(itemSeparator) => handlePrefix(values.map(valueMapper(_).valueString).mkString(itemSeparator))

        /*
        via: http://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineBinding

        > If itemSeparator is specified, add prefix and the join the array into a single string with itemSeparator
        > separating the items. Otherwise first add prefix, then recursively process individual elements.

        Not 100% sure if this is conformant, as we are only recurse the head into `processValue` here... However the
        conformance test "Test InlineJavascriptRequirement with multiple expressions in the same tool" is happy with the
        behavior.
         */
        // InstantiatedCommand elements are generated here when there exists an optionalValueFrom
        case None if commandLineBinding.optionalValueFrom.isDefined => values.toList match {
          case head :: tail => processValue(head) ++ tail.map(valueMapper(_).valueString)
          case Nil => prefixAsList
        }

        // When there is no optionalValueFrom the InstantiatedCommand elements for the womValue are appended elsewhere
        case _ => prefixAsList
      }
      case _: WomObjectLike => prefixAsList
      case WomEnumerationValue(_, value) => handlePrefix(value)
      case WomCoproductValue(_, value) => processValue(value)
      case w => throw new RuntimeException(s"Unhandled CwlExpressionCommandPart value '$w' of type ${w.womType.stableName}")
    }

    evaluatedWomValue map { v => processValue(v) map applyShellQuote map (InstantiatedCommand(_)) } toValidated
  }
}

case class InputCommandLineBindingCommandPart(commandLineBinding: InputCommandLineBinding, associatedValue: WomValue)
                                             (hasShellCommandRequirement: Boolean, expressionLib: ExpressionLib)
  extends CommandLineBindingCommandPart(commandLineBinding)(hasShellCommandRequirement, expressionLib) {
  override lazy val boundValue = Option(associatedValue)
}

case class ArgumentCommandLineBindingCommandPart(commandLineBinding: ArgumentCommandLineBinding)
                                                (hasShellCommandRequirement: Boolean, expressionLib: ExpressionLib)
  extends CommandLineBindingCommandPart(commandLineBinding)(hasShellCommandRequirement, expressionLib) {
  override lazy val boundValue = None
}
