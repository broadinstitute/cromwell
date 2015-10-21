package cromwell.binding.expression

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.WdlExpression._
import cromwell.binding.types._
import cromwell.binding.{WdlExpressionException, WdlNamespace}
import cromwell.parser.WdlParser.{Ast, AstNode, Terminal}
import cromwell.util.TryUtil

import scala.util.{Failure, Success, Try}

case class TypeEvaluator(override val lookup: String => WdlType, override val functions: WdlFunctions[WdlType]) extends Evaluator {
  override type T = WdlType

  override def evaluate(ast: AstNode): Try[WdlType] = ast match {
    case t: Terminal if t.getTerminalStr == "identifier" => Try(lookup(t.getSourceString))
    case t: Terminal if t.getTerminalStr == "integer" => Success(WdlIntegerType)
    case t: Terminal if t.getTerminalStr == "float" => Success(WdlFloatType)
    case t: Terminal if t.getTerminalStr == "boolean" => Success(WdlBooleanType)
    case t: Terminal if t.getTerminalStr == "string" => Success(WdlStringType)
    case a: Ast if a.isBinaryOperator =>
      val lhs = evaluate(a.getAttribute("lhs"))
      val rhs = evaluate(a.getAttribute("rhs"))
      a.getName match {
        case "Add" => for(l <- lhs; r <- rhs) yield l.add(r).get
        case "Subtract" => for(l <- lhs; r <- rhs) yield l.subtract(r).get
        case "Multiply" => for(l <- lhs; r <- rhs) yield l.multiply(r).get
        case "Divide" => for(l <- lhs; r <- rhs) yield l.divide(r).get
        case "Remainder" => for(l <- lhs; r <- rhs) yield l.mod(r).get
        case "Equals" => for(l <- lhs; r <- rhs) yield l.equals(r).get
        case "NotEquals" => for(l <- lhs; r <- rhs) yield l.notEquals(r).get
        case "LessThan" => for(l <- lhs; r <- rhs) yield l.lessThan(r).get
        case "LessThanOrEqual" => for(l <- lhs; r <- rhs) yield l.lessThanOrEqual(r).get
        case "GreaterThan" => for(l <- lhs; r <- rhs) yield l.greaterThan(r).get
        case "GreaterThanOrEqual" => for(l <- lhs; r <- rhs) yield l.greaterThanOrEqual(r).get
        case "LogicalOr" => for(l <- lhs; r <- rhs) yield l.or(r).get
        case "LogicalAnd" => for(l <- lhs; r <- rhs) yield l.and(r).get
        case _ => Failure(new WdlExpressionException(s"Invalid operator: ${a.getName}"))
      }
    case a: Ast if a.isUnaryOperator =>
      val expression = evaluate(a.getAttribute("expression"))
      a.getName match {
        case "LogicalNot" => for(e <- expression) yield e.not.get
        case "UnaryPlus" => for(e <- expression) yield e.unaryPlus.get
        case "UnaryNegation" => for(e <- expression) yield e.unaryMinus.get
        case _ => Failure(new WdlExpressionException(s"Invalid operator: ${a.getName}"))
      }
    case a: Ast if a.isArrayLiteral =>
      val evaluatedElements = a.getAttribute("values").astListAsVector map evaluate
      for {
        flattened <- TryUtil.flatten(evaluatedElements)
        subtype <- WdlType.homogeneousTypeFromTypes(flattened)
      } yield WdlArrayType(subtype)
    case a: Ast if a.isMapLiteral =>
      val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
        val key = evaluate(kv.asInstanceOf[Ast].getAttribute("key"))
        val value = evaluate(kv.asInstanceOf[Ast].getAttribute("value"))
        key -> value
      }

      val flattenedTries = evaluatedMap flatMap { case (k,v) => Seq(k,v) }
      flattenedTries partition {_.isSuccess} match {
        case (_, failures) if failures.nonEmpty =>
          val message = failures.collect { case f: Failure[_] => f.exception.getMessage }.mkString("\n")
          Failure(new WdlExpressionException(s"Could not evaluate expression:\n$message"))
        case (successes, _) =>
          for {
            keyType <- WdlType.homogeneousTypeFromTypes(evaluatedMap map { case (k, v) => k.get})
            valueType <- WdlType.homogeneousTypeFromTypes(evaluatedMap map { case (k, v) => v.get})
          } yield WdlMapType(keyType, valueType)
      }
    case a: Ast if a.isMemberAccess =>
      a.getAttribute("rhs") match {
        case rhs: Terminal if rhs.getTerminalStr == "identifier" =>
          evaluate(a.getAttribute("lhs")).flatMap {
            case o: WdlCallOutputsObjectType =>
              o.call.task.outputs.find(_.name == rhs.getSourceString) match {
                case Some(taskOutput) => evaluate(taskOutput.expression.ast)
                case None => Failure(new WdlExpressionException(s"Could not find key ${rhs.getSourceString}"))
              }
            case ns: WdlNamespace => Success(lookup(ns.importedAs.map {n => s"$n.${rhs.getSourceString}"}.getOrElse(rhs.getSourceString)))
            case _ => Failure(new WdlExpressionException("Left-hand side of expression must be a WdlObject or Namespace"))
          }
        case _ => Failure(new WdlExpressionException("Right-hand side of expression must be identifier"))
      }
    case a: Ast if a.isArrayOrMapLookup =>
      (evaluate(a.getAttribute("lhs")), evaluate(a.getAttribute("rhs"))) match {
        case (Success(a: WdlArrayType), Success(WdlIntegerType)) => Success(a.memberType)
        case (Success(m: WdlMapType), Success(v: WdlType)) => Success(m.valueType)
        case (Failure(ex), _) => Failure(ex)
        case (_, Failure(ex)) => Failure(ex)
        case (_, _) => Failure(new WdlExpressionException(s"Can't index ${a.toPrettyString}"))
      }
    case a: Ast if a.isFunctionCall =>
      val name = a.getAttribute("name").sourceString
      val params = a.params map evaluate
      functions.getFunction(name)(params)
  }
}

