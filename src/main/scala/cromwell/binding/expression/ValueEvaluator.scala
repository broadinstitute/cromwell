package cromwell.binding.expression

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.WdlExpression._
import cromwell.binding.types._
import cromwell.binding.values.{WdlValue, _}
import cromwell.binding.{WdlExpressionException, WdlNamespace}
import cromwell.parser.WdlParser.{Ast, AstNode, Terminal}

import scala.util.{Failure, Success, Try}

case class ValueEvaluator(override val lookup: String => WdlValue, override val functions: WdlFunctions[WdlValue]) extends Evaluator {
  override type T = WdlValue

  private def replaceInterpolationTag(string: String, tag: String, lookup: ScopedLookupFunction) =
    string.replace(tag, lookup(tag.substring(2, tag.length - 1)).valueString)

  private def interpolate(str: String, lookup: ScopedLookupFunction): String =
    "\\$\\{([a-zA-Z]([a-zA-Z0-9_])*)\\}".r.findAllIn(str).foldLeft(str) {replaceInterpolationTag(_, _, lookup)}

  override def evaluate(ast: AstNode): Try[WdlValue] = {
    ast match {
      case t: Terminal if t.getTerminalStr == "identifier" => Try(lookup(t.getSourceString))
      case t: Terminal if t.getTerminalStr == "integer" => Success(WdlInteger(t.getSourceString.toInt))
      case t: Terminal if t.getTerminalStr == "float" => Success(WdlFloat(t.getSourceString.toDouble))
      case t: Terminal if t.getTerminalStr == "boolean" => Success(WdlBoolean(t.getSourceString == "true"))
      case t: Terminal if t.getTerminalStr == "string" =>
        Try(WdlString(interpolate(t.getSourceString, lookup)))
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
        val evaluatedElements = a.getAttribute("values").astListAsVector.map(evaluate)
        evaluatedElements.partition {_.isSuccess} match {
          case (_, failures) if failures.nonEmpty =>
            val message = failures.collect {case f: Failure[_] => f.exception.getMessage}.mkString("\n")
            Failure(new WdlExpressionException(s"Could not evaluate expression:\n$message"))
          case (successes, _) =>
            for (subtype <- WdlType.homogeneousTypeFromValues(successes.map(_.get)))
              yield WdlArray(WdlArrayType(subtype), successes.map(_.get))
        }
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
            WdlMapType(WdlAnyType, WdlAnyType).coerceRawValue(evaluatedMap.map({ case (k, v) => k.get -> v.get }).toMap)
        }
      case a: Ast if a.isMemberAccess =>
        a.getAttribute("rhs") match {
          case rhs:Terminal if rhs.getTerminalStr == "identifier" =>
            evaluate(a.getAttribute("lhs")).flatMap {
              case o: WdlObjectLike =>
                o.value.get(rhs.getSourceString) match {
                  case Some(v:WdlValue) => Success(v)
                  case None => Failure(new WdlExpressionException(s"Could not find key ${rhs.getSourceString}"))
                }
              case a: WdlArray if a.wdlType == WdlArrayType(WdlObjectType) =>
                /**
                 * This case is for slicing an Array[Object], used mainly for scatter-gather.
                 * For example, if 'call foo' was in a scatter block, foo's outputs (e.g. Int x)
                 * would be an Array[Int].  If a downstream call has an input expression "foo.x",
                 * then 'foo' would evaluate to an Array[Objects] and foo.x would result in an
                 * Array[Int]
                 */
                Success(a map {_.asInstanceOf[WdlObject].value.get(rhs.sourceString).get})
              case ns: WdlNamespace => Success(lookup(ns.importedAs.map {n => s"$n.${rhs.getSourceString}"}.getOrElse(rhs.getSourceString)))
              case _ => Failure(new WdlExpressionException("Left-hand side of expression must be a WdlObject or Namespace"))
            }
          case _ => Failure(new WdlExpressionException("Right-hand side of expression must be identifier"))
        }
      case a: Ast if a.isArrayOrMapLookup =>
        val index = evaluate(a.getAttribute("rhs"))
        val mapOrArray = evaluate(a.getAttribute("lhs"))
        (mapOrArray, index) match {
          case (Success(a: WdlArray), Success(i: WdlInteger)) =>
            Try(a.value(i.value)) match {
              case s:Success[WdlValue] => s
              case Failure(ex) => Failure(new WdlExpressionException(s"Failed to find index $index on array:\n\n$mapOrArray\n\n${ex.getMessage}"))
            }
          case (Success(m: WdlMap), Success(v: WdlValue)) =>
            m.value.get(v) match {
              case Some(value) => Success(value)
              case _ => Failure(new WdlExpressionException(s"Failed to find a key '$index' on a map:\n\n$mapOrArray"))
            }
          case (Failure(ex), _) => Failure(ex)
          case (_, Failure(ex)) => Failure(ex)
          case (_, _) => Failure(new WdlExpressionException(s"Can't index $mapOrArray with index $index"))
        }
      case a: Ast if a.isFunctionCall =>
        val name = a.getAttribute("name").sourceString
        val params = a.params map evaluate
        functions.getFunction(name)(params)
    }
  }
}

