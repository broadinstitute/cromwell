package wdl4s.expression

import lenthall.util.TryUtil
import wdl4s.AstTools.EnhancedAstNode
import wdl4s.WdlExpression._
import wdl4s.types._
import wdl4s.values.{WdlValue, _}
import wdl4s._
import wdl4s.parser.WdlParser.{Ast, AstNode, Terminal}
import wdl4s.TernaryIf

import scala.util.{Failure, Success, Try}

case class ValueEvaluator(override val lookup: String => WdlValue, override val functions: WdlFunctions[WdlValue]) extends Evaluator {
  override type T = WdlValue

  private val InterpolationTagPattern = "\\$\\{\\s*([^\\}]*)\\s*\\}".r

  private def replaceInterpolationTag(string: Try[WdlString], tag: String): Try[WdlString] = {
    val expr = WdlExpression.fromString(tag.substring(2, tag.length - 1))
    (expr.evaluate(lookup, functions), string) match {
      case (Success(value), Success(str)) =>
        value match {
          case s: WdlString if InterpolationTagPattern.anchored.findFirstIn(s.valueString).isDefined =>
            replaceInterpolationTag(Success(WdlString(str.value.replace(tag, s.valueString))), s.valueString)
          case v => Success(WdlString(str.value.replace(tag, v.valueString)))
        }
      case (Failure(ex), _) => Failure(ex)
      case (_, Failure(ex)) => Failure(ex)
    }
  }

  private def interpolate(str: String): Try[WdlString] = {
    InterpolationTagPattern.findAllIn(str).foldLeft(Try(WdlString(str)))(replaceInterpolationTag)
  }

  private def interpolate(value: WdlValue): Try[WdlValue] = {
    value match {
      case s: WdlString => interpolate(s.valueString)
      case _ => Try(value)
    }
  }

  override def evaluate(ast: AstNode): Try[WdlValue] = {
    ast match {
      case t: Terminal if t.getTerminalStr == "identifier" => Try(lookup(t.getSourceString)).flatMap(interpolate)
      case t: Terminal if t.getTerminalStr == "integer" => Success(WdlInteger(t.getSourceString.toInt))
      case t: Terminal if t.getTerminalStr == "float" => Success(WdlFloat(t.getSourceString.toDouble))
      case t: Terminal if t.getTerminalStr == "boolean" => Success(WdlBoolean(t.getSourceString == "true"))
      case t: Terminal if t.getTerminalStr == "string" => interpolate(t.getSourceString)
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
      case TernaryIf(condition, ifTrue, ifFalse) =>
        evaluate(condition) flatMap {
          case WdlBoolean(true) => evaluate(ifTrue)
          case WdlBoolean(false) => evaluate(ifFalse)
          case other => Failure(new WdlExpressionException("'if' expression must be given a boolean argument but got: " + other.toWdlString))
        }
      case a: Ast if a.isArrayLiteral =>
        val evaluatedElements = a.getAttribute("values").astListAsVector map evaluate
        for {
          elements <- TryUtil.sequence(evaluatedElements)
          subtype = WdlType.homogeneousTypeFromValues(elements)
          isEmpty = elements.isEmpty
        } yield WdlArray(WdlArrayType(subtype), elements.map(subtype.coerceRawValue(_).get))
      case a: Ast if a.isTupleLiteral =>
        val unevaluatedElements = a.getAttribute("values").astListAsVector
        if (unevaluatedElements.size == 1) {
          evaluate(unevaluatedElements.head)
        } else if (unevaluatedElements.size == 2) {
          for {
            left <- evaluate(unevaluatedElements.head)
            right <- evaluate(unevaluatedElements(1))
          } yield WdlPair(left, right)
        } else {
          Failure(new WdlExpressionException(s"WDL does not currently support tuples with n > 2: ${a.toPrettyString}"))
        }
      case a: Ast if a.isMapLiteral =>
        val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
          val key = evaluate(kv.asInstanceOf[Ast].getAttribute("key"))
          val value = evaluate(kv.asInstanceOf[Ast].getAttribute("value"))
          key -> value
        }
        TryUtil.sequence(evaluatedMap map { tuple => TryUtil.sequenceTuple(tuple) }) flatMap { pairs =>
          WdlMapType(WdlAnyType, WdlAnyType).coerceRawValue(pairs.toMap)
        }
      case a: Ast if a.isMemberAccess =>
        a.getAttribute("rhs") match {
          case rhs:Terminal if rhs.getTerminalStr == "identifier" =>
            evaluate(a.getAttribute("lhs")).flatMap {
              case o: WdlObjectLike =>
                o.value.get(rhs.getSourceString) match {
                  case Some(v:WdlValue) => Success(v)
                  case None =>
                    o match {
                        // o is a CallOutputsObject which means we failed to find an output value for rhs
                        // Give a specific error message based on the type of Callable
                      case callOutputObject: WdlCallOutputsObject => 
                        callOutputObject.call match {
                          case workflowCall: WorkflowCall => 
                            Failure(new WdlExpressionException(
                              s"""${rhs.getSourceString} is not declared as an output of the sub workflow ${workflowCall.calledWorkflow.fullyQualifiedName}.
                                 |If you want to use workflow ${workflowCall.calledWorkflow.fullyQualifiedName} as a sub workflow, make sure that its output section is up to date with the latest syntax.
                                 |See the WDL specification for how to write outputs: https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs""".stripMargin
                            ))
                          case taskCall: TaskCall => 
                            Failure(new WdlExpressionException(
                              s"""${rhs.getSourceString} is not declared as an output of the task ${taskCall.task.fullyQualifiedName}.
                                 |Make sure to declare it as an output to be able to use it in the workflow.""".stripMargin
                            ))
                          case unknownCall => 
                            Failure(new WdlExpressionException(
                              s"Could not find key ${rhs.getSourceString} in Call ${unknownCall.fullyQualifiedName} of unknown type."
                            ))
                        }
                      case _ => Failure(new WdlExpressionException(s"Could not find key ${rhs.getSourceString} in WdlObject"))
                    }
                }
              case array: WdlArray if array.wdlType == WdlArrayType(WdlObjectType) =>
                /**
                 * This case is for slicing an Array[Object], used mainly for scatter-gather.
                 * For example, if 'call foo' was in a scatter block, foo's outputs (e.g. Int x)
                 * would be an Array[Int].  If a downstream call has an input expression "foo.x",
                 * then 'foo' would evaluate to an Array[Objects] and foo.x would result in an
                 * Array[Int]
                 */
                Success(array map {_.asInstanceOf[WdlObject].value.get(rhs.sourceString).get})
              case p: WdlPair =>
                val identifier = rhs.getSourceString
                if (identifier.equals("left")) Success(p.left)
                else if (identifier.equals("right")) Success(p.right)
                else Failure(new WdlExpressionException("A pair only has the members: 'left' and 'right'"))
              case ns: WdlNamespace => Success(lookup(ns.importedAs.map{ n => s"$n.${rhs.getSourceString}" }.getOrElse(rhs.getSourceString)))
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

