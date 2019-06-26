package wdl.draft2.model.expression

import common.util.TryUtil
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.model.WdlExpression._
import wdl.draft2.model._
import wdl.draft2.model.values.WdlCallOutputsObject
import wdl.draft2.parser.WdlParser.{Ast, AstNode, Terminal}
import wom.WomExpressionException
import wom.types.{WomArrayType, WomMapType, WomObjectType, WomType}
import wom.values._

import scala.util.{Failure, Success, Try}

object ValueEvaluator {
  val InterpolationTagPattern = "\\$\\{\\s*([^\\}]*)\\s*\\}".r
}

case class ValueEvaluator(override val lookup: String => WomValue, override val functions: WdlFunctions[WomValue]) extends Evaluator {
  override type T = WomValue

  private val InterpolationTagPattern = "\\$\\{\\s*([^\\}]*)\\s*\\}".r


  private def interpolate(strToProcess: String, resultSoFar: Try[WomString] = Success(WomString(""))): Try[WomString] = {

    def evaluateTag(tag: String): Try[WomString] = {
      val expr = WdlExpression.fromString(tag.substring(2, tag.length - 1))
      expr.evaluate(lookup, functions).map(result => WomString(result.valueString))
    }

    InterpolationTagPattern.findFirstIn(strToProcess) match {
      case Some(tag) =>
        val prefix = strToProcess.substring(0, strToProcess.indexOf(tag))
        val suffix = strToProcess.substring(strToProcess.indexOf(tag) + tag.length, strToProcess.length)

        val newResultSoFar = for {
          rsf <- resultSoFar
          interpolatedTag <- evaluateTag(tag)
        } yield WomString(rsf.value + prefix + interpolatedTag.value)

        interpolate(suffix, newResultSoFar)

      case None => resultSoFar.map(ws => WomString(ws.value + strToProcess))
    }

  }

  override def evaluate(ast: AstNode): Try[WomValue] = {
    ast match {
      case t: Terminal if t.getTerminalStr == "identifier" => Try(lookup(t.getSourceString))
      case t: Terminal if t.getTerminalStr == "integer" => Success(WomInteger(t.getSourceString.toInt))
      case t: Terminal if t.getTerminalStr == "float" => Success(WomFloat(t.getSourceString.toDouble))
      case t: Terminal if t.getTerminalStr == "boolean" => Success(WomBoolean(t.getSourceString == "true"))
      case t: Terminal if t.getTerminalStr == "string" => interpolate(t.getSourceString)
      case a: Ast if a.isBinaryOperator =>
        val lhs = evaluate(a.getAttribute("lhs"))
        val rhs = () => evaluate(a.getAttribute("rhs"))
        a.getName match {
          case "Add" => for(l <- lhs; r <- rhs()) yield l.add(r).get
          case "Subtract" => for(l <- lhs; r <- rhs()) yield l.subtract(r).get
          case "Multiply" => for(l <- lhs; r <- rhs()) yield l.multiply(r).get
          case "Divide" => for(l <- lhs; r <- rhs()) yield l.divide(r).get
          case "Remainder" => for(l <- lhs; r <- rhs()) yield l.mod(r).get
          case "Equals" => for(l <- lhs; r <- rhs()) yield l.equals(r).get
          case "NotEquals" => for(l <- lhs; r <- rhs()) yield l.notEquals(r).get
          case "LessThan" => for(l <- lhs; r <- rhs()) yield l.lessThan(r).get
          case "LessThanOrEqual" => for(l <- lhs; r <- rhs()) yield l.lessThanOrEqual(r).get
          case "GreaterThan" => for(l <- lhs; r <- rhs()) yield l.greaterThan(r).get
          case "GreaterThanOrEqual" => for(l <- lhs; r <- rhs()) yield l.greaterThanOrEqual(r).get
          case "LogicalOr" => lhs flatMap {
            case WomBoolean(true) => Success(WomBoolean(true))
            case b => rhs() flatMap b.or
          }
          case "LogicalAnd" => lhs flatMap {
            case WomBoolean(false) => Success(WomBoolean(false))
            case b => rhs() flatMap b.and
          }
          case _ => Failure(new WomExpressionException(s"Invalid operator: ${a.getName}"))
        }
      case a: Ast if a.isUnaryOperator =>
        val expression = evaluate(a.getAttribute("expression"))
        a.getName match {
          case "LogicalNot" => for(e <- expression) yield e.not.get
          case "UnaryPlus" => for(e <- expression) yield e.unaryPlus.get
          case "UnaryNegation" => for(e <- expression) yield e.unaryMinus.get
          case _ => Failure(new WomExpressionException(s"Invalid operator: ${a.getName}"))
        }
      case TernaryIf(condition, ifTrue, ifFalse) =>
        evaluate(condition) flatMap {
          case WomBoolean(true) => evaluate(ifTrue)
          case WomBoolean(false) => evaluate(ifFalse)
          case other => Failure(new WomExpressionException("'if' expression must be given a boolean argument but got: " + other.toWomString))
        }
      case a: Ast if a.isArrayLiteral =>
        val evaluatedElements = a.getAttribute("values").astListAsVector map evaluate
        for {
          elements <- TryUtil.sequence(evaluatedElements)
          subtype = WomType.homogeneousTypeFromValues(elements)
          isEmpty = elements.isEmpty
        } yield WomArray(WomArrayType(subtype), elements.map(subtype.coerceRawValue(_).get))
      case a: Ast if a.isTupleLiteral =>
        val unevaluatedElements = a.getAttribute("values").astListAsVector
        if (unevaluatedElements.size == 1) {
          evaluate(unevaluatedElements.head)
        } else if (unevaluatedElements.size == 2) {
          for {
            left <- evaluate(unevaluatedElements.head)
            right <- evaluate(unevaluatedElements(1))
          } yield WomPair(left, right)
        } else {
          Failure(new WomExpressionException(s"WDL does not currently support tuples with n > 2: ${a.toPrettyString}"))
        }
      case a: Ast if a.isMapLiteral =>
        val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
          val key = evaluate(kv.asInstanceOf[Ast].getAttribute("key"))
          val value = evaluate(kv.asInstanceOf[Ast].getAttribute("value"))
          key -> value
        }
        TryUtil.sequence(evaluatedMap map { tuple => TryUtil.sequenceTuple(tuple) }) flatMap { pairs =>
          WomMapType(
            WomType.homogeneousTypeFromValues(pairs.map(_._1)),
            WomType.homogeneousTypeFromValues(pairs.map(_._2))
          ).coerceRawValue(pairs.toMap)
        }
      case a: Ast if a.isObjectLiteral =>
        val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
          val key = kv.asInstanceOf[Ast].getAttribute("key").sourceString
          val value = evaluate(kv.asInstanceOf[Ast].getAttribute("value"))
          key -> value
        }
        TryUtil.sequence(evaluatedMap map { tuple => tuple._2.map((tuple._1, _)) }) map { pairs =>
          WomObject(pairs.toMap)
        }
      case a: Ast if a.isMemberAccess =>
        a.getAttribute("rhs") match {
          case rhs:Terminal if rhs.getTerminalStr == "identifier" =>
            val memberAccessAsString = s"${a.getAttribute("lhs").sourceString}.${a.getAttribute("rhs").sourceString}"
            Try(lookup(memberAccessAsString)).recoverWith {
              case _ =>
                evaluate(a.getAttribute("lhs")).flatMap {
                  case o: WomObjectLike =>
                    o.values.get(rhs.getSourceString) match {
                      case Some(v:WomValue) => Success(v)
                      case None =>
                        o match {
                          // o is a CallOutputsObject which means we failed to find an output value for rhs
                          // Give a specific error message based on the type of Callable
                          case callOutputObject: WdlCallOutputsObject =>
                            callOutputObject.call match {
                              case workflowCall: WdlWorkflowCall =>
                                Failure(new WomExpressionException(
                                  s"""${rhs.getSourceString} is not declared as an output of the sub workflow ${workflowCall.calledWorkflow.fullyQualifiedName}.
                                     |If you want to use workflow ${workflowCall.calledWorkflow.fullyQualifiedName} as a sub workflow, make sure that its output section is up to date with the latest syntax.
                                     |See the WDL specification for how to write outputs: https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs""".stripMargin
                                ))
                              case taskCall: WdlTaskCall =>
                                Failure(new WomExpressionException(
                                  s"""${rhs.getSourceString} is not declared as an output of the task ${taskCall.task.fullyQualifiedName}.
                                     |Make sure to declare it as an output to be able to use it in the workflow.""".stripMargin
                                ))
                              case unknownCall =>
                                Failure(new WomExpressionException(
                                  s"Could not find key ${rhs.getSourceString} in Call ${unknownCall.fullyQualifiedName} of unknown type."
                                ))
                            }
                          case _ => Failure(new WomExpressionException(s"Could not find key ${rhs.getSourceString} in WdlObject"))
                        }
                    }
                  case array: WomArray if array.womType == WomArrayType(WomObjectType) =>
                    /*
                     * This case is for slicing an Array[Object], used mainly for scatter-gather.
                     * For example, if 'call foo' was in a scatter block, foo's outputs (e.g. Int x)
                     * would be an Array[Int].  If a downstream call has an input expression "foo.x",
                     * then 'foo' would evaluate to an Array[Objects] and foo.x would result in an
                     * Array[Int]
                     */
                    Success(array map {_.asInstanceOf[WomObject].values.get(rhs.sourceString).get})
                  case p: WomPair =>
                    val identifier = rhs.getSourceString
                    if (identifier.equals("left")) Success(p.left)
                    else if (identifier.equals("right")) Success(p.right)
                    else Failure(new WomExpressionException("A pair only has the members: 'left' and 'right'"))
                  case ns: WdlNamespace => Success(lookup(ns.importedAs.map{ n => s"$n.${rhs.getSourceString}" }.getOrElse(rhs.getSourceString)))
                  case _ => Failure(new WomExpressionException("Left-hand side of expression must be a WdlObject or Namespace"))
                }
            }
          case _ => Failure(new WomExpressionException("Right-hand side of expression must be identifier"))
        }
      case a: Ast if a.isArrayOrMapLookup =>
        val index = evaluate(a.getAttribute("rhs"))
        val mapOrArray = evaluate(a.getAttribute("lhs"))
        (mapOrArray, index) match {
          case (Success(a: WomArray), Success(i: WomInteger)) =>
            Try(a.value(i.value)) match {
              case s:Success[WomValue] => s
              case Failure(ex) => Failure(new WomExpressionException(s"Failed to find index $index on array:\n\n$mapOrArray\n\n${ex.getMessage}"))
            }
          case (Success(m: WomMap), Success(v: WomValue)) =>
            m.value.get(v) match {
              case Some(value) => Success(value)
              case _ => Failure(new WomExpressionException(s"Failed to find a key '$index' on a map:\n\n$mapOrArray"))
            }
          case (Failure(ex), _) => Failure(ex)
          case (_, Failure(ex)) => Failure(ex)
          case (_, _) => Failure(new WomExpressionException(s"Can't index $mapOrArray with index $index"))
        }
      case a: Ast if a.isFunctionCall =>
        val name = a.getAttribute("name").sourceString
        val params = a.params map evaluate
        functions.getFunction(name)(params)
    }
  }
}

