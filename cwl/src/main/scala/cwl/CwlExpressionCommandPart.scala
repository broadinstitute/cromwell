package cwl

import cats.data.Validated.Invalid
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

case class CwlExpressionCommandPart(expr: Expression)(expressionLib: ExpressionLib) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
    val stringKeyMap = inputsMap map {
      case (LocalName(localName), value) => localName -> value
    }
    val parameterContext = ParameterContext(
        inputs = stringKeyMap, runtimeOption = Option(runtimeEnvironment)
      )
    ExpressionEvaluator.eval(expr, parameterContext, expressionLib) map { womValue =>
      List(InstantiatedCommand(womValue.valueString))
    }
  }
}

/**
  * Generates command parts from a CWL Binding
  */
abstract class CommandLineBindingCommandPart(commandLineBinding: CommandLineBinding)(expressionLib: ExpressionLib) extends CommandPart {

  
  private lazy val prefixAsString = commandLineBinding.prefix.getOrElse("")
  private lazy val prefixAsList = commandLineBinding.prefix.toList
  private lazy val separate = commandLineBinding.effectiveSeparate

  def handlePrefix(value: String) = {
    if (separate) prefixAsList :+ value else List(s"$prefixAsString$value")
  }

  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
      val stringInputsMap = inputsMap map {
        case (LocalName(localName), value) => localName -> valueMapper(value)
      }
      val parameterContext = ParameterContext(inputs = stringInputsMap, runtimeOption = Option(runtimeEnvironment))

    val evaluatedValueFrom = commandLineBinding.optionalValueFrom map {
      case StringOrExpression.Expression(expression) => ExpressionEvaluator.eval(expression, parameterContext, expressionLib) map valueMapper
      case StringOrExpression.String(string) => WomString(string).validNel
    }

    val evaluatedWomValue: Checked[WomValue] = evaluatedValueFrom.orElse(boundValue.map(_.validNel)) match {
      case Some(womValue) => womValue.map(valueMapper).toEither
      case None => "Command line binding has no valueFrom field and no bound value".invalidNelCheck
    }
    
    def applyShellQuote(womValue: WomValue): Checked[WomValue] = commandLineBinding.shellQuote match {
      case Some(false) | None => womValue.validNelCheck
      case _ => "Shell quote is unsupported yet".invalidNelCheck
    }
    
    def processValue(womValue: WomValue): List[String] = womValue match {
      case WomOptionalValue(_, Some(womValue)) => processValue(womValue)
      case _: WomString | _: WomInteger | _: WomFile => handlePrefix(womValue.valueString)
      // For boolean values, use the value of the boolean to choose whether to print the prefix or not
      case WomBoolean(false) => List.empty
      case WomBoolean(true) => prefixAsList
      case WomArray(_, values) => commandLineBinding.itemSeparator match {
        case Some(itemSeparator) => handlePrefix(values.map(_.valueString).mkString(itemSeparator))
        case None if commandLineBinding.optionalValueFrom.isDefined => values.toList.flatMap(processValue)
        case _ => prefixAsList
      }
      case _: WomObjectLike => prefixAsList
      case _ => List.empty
    }

    // Can't flatMap Either in 2.11.
    evaluatedWomValue match {
      case Right(womValue) => applyShellQuote(womValue).map(processValue).map(_.map(InstantiatedCommand(_))).toValidated
      case Left(e) => Invalid(e)
    }
  }
  
  def boundValue: Option[WomValue]
}

case class InputCommandLineBindingCommandPart(commandLineBinding: InputCommandLineBinding, associatedValue: WomValue)(expressionLib: ExpressionLib) extends CommandLineBindingCommandPart(commandLineBinding)(expressionLib) {
  override lazy val boundValue = Option(associatedValue)
}

case class ArgumentCommandLineBindingCommandPart(commandLineBinding: ArgumentCommandLineBinding)(expressionLib: ExpressionLib) extends CommandLineBindingCommandPart(commandLineBinding)(expressionLib) {
  override lazy val boundValue = None
}
