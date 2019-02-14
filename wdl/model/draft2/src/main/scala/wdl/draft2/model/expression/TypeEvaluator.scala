package wdl.draft2.model.expression

import common.util.TryUtil
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.model.WdlExpression._
import wdl.draft2.model.types.WdlCallOutputsObjectType
import wdl.draft2.model.{DeclarationInterface, Scope, TernaryIf, WdlNamespace}
import wdl.draft2.parser.WdlParser.{Ast, AstNode, Terminal}
import wom.WomExpressionException
import wom.types._

import scala.util.{Failure, Success, Try}

case class TypeEvaluator(override val lookup: String => WomType, override val functions: WdlFunctions[WomType], from: Option[Scope] = None) extends Evaluator {
  override type T = WomType

  override def evaluate(ast: AstNode): Try[WomType] = ast match {
    case null => Failure(new Exception("Cannot evaluate the type of an empty expression"))
    case t: Terminal if t.getTerminalStr == "identifier" => Try(lookup(t.getSourceString))
    case t: Terminal if t.getTerminalStr == "integer" => Success(WomIntegerType)
    case t: Terminal if t.getTerminalStr == "float" => Success(WomFloatType)
    case t: Terminal if t.getTerminalStr == "boolean" => Success(WomBooleanType)
    case t: Terminal if t.getTerminalStr == "string" => Success(WomStringType)
    case a: Ast if a.isBinaryOperator =>
      val lhs = evaluate(a.getAttribute("lhs"))
      val rhs = evaluate(a.getAttribute("rhs"))
      a.getName match {
        case "Add" => for(l <- lhs; r <- rhs) yield l.add(r).get
        case "Subtract" => for(l <- lhs; r <- rhs) yield l.subtract(r).get
        case "Multiply" => for(l <- lhs; r <- rhs) yield l.multiply(r).get
        case "Divide" => for(l <- lhs; r <- rhs) yield l.divide(r).get
        case "Remainder" => for(l <- lhs; r <- rhs) yield l.mod(r).get
        case "Equals" => for(l <- lhs; r <- rhs) yield l.equalsType(r).get
        case "NotEquals" => for(l <- lhs; r <- rhs) yield l.notEquals(r).get
        case "LessThan" => for(l <- lhs; r <- rhs) yield l.lessThan(r).get
        case "LessThanOrEqual" => for(l <- lhs; r <- rhs) yield l.lessThanOrEqual(r).get
        case "GreaterThan" => for(l <- lhs; r <- rhs) yield l.greaterThan(r).get
        case "GreaterThanOrEqual" => for(l <- lhs; r <- rhs) yield l.greaterThanOrEqual(r).get
        case "LogicalOr" => for(l <- lhs; r <- rhs) yield l.or(r).get
        case "LogicalAnd" => for(l <- lhs; r <- rhs) yield l.and(r).get
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
        case WomBooleanType => for {
          ifTrueType <- evaluate(ifTrue)
          ifFalseType <- evaluate(ifFalse)
        } yield WomType.lowestCommonSubtype(Seq(ifTrueType, ifFalseType))
        case _ => Failure(new WomExpressionException("The condition of a ternary 'if' must be a Boolean."))
      }
    case a: Ast if a.isArrayLiteral =>
      val evaluatedElements = a.getAttribute("values").astListAsVector map evaluate
      for {
        elements <- TryUtil.sequence(evaluatedElements)
        subtype = WomType.homogeneousTypeFromTypes(elements)
      } yield WomArrayType(subtype)
    case o: Ast if o.isObjectLiteral => Success(WomObjectType)
    case a: Ast if a.isTupleLiteral => tupleAstToWdlType(a)
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
          Failure(new WomExpressionException(s"Could not evaluate expression:\n$message"))
        case good @ _ =>
          val keyType = WomType.homogeneousTypeFromTypes(evaluatedMap map { case (k, _) => k.get} )
          val valueType = WomType.homogeneousTypeFromTypes(evaluatedMap map { case (_, v) => v.get} )
          Success(WomMapType(keyType, valueType))
      }
    case a: Ast if a.isMemberAccess =>
      a.getAttribute("rhs") match {
        case rhs: Terminal if rhs.getTerminalStr == "identifier" =>
          evaluate(a.getAttribute("lhs")).flatMap {
            case o: WdlCallOutputsObjectType =>
              o.call.outputs.find(_.unqualifiedName == rhs.getSourceString) match {
                case Some(taskOutput) =>
                  val t = taskOutput.womType
                  val relative = from match {
                    case None => t
                    case Some(scope) => DeclarationInterface.relativeWdlType(scope, taskOutput, t)
                  }
                  Success(relative)
                case None => Failure(new WomExpressionException(s"Could not find key ${rhs.getSourceString}"))
              }
            case WomPairType(leftType, rightType) =>
              rhs.sourceString match {
                case "left" => Success(leftType)
                case "right" => Success(rightType)
              }
            case WomObjectType => Success(WomAnyType)
            case ns: WdlNamespace => Success(lookup(ns.importedAs.map{ n => s"$n.${rhs.getSourceString}" }.getOrElse(rhs.getSourceString)))
            case _ => Failure(new WomExpressionException("Left-hand side of expression must be a WdlObject or Namespace"))
          } recoverWith {
            case _ => Try(lookup(a.getAttribute("lhs").sourceString + "." + rhs.sourceString))
          }
        case _ => Failure(new WomExpressionException("Right-hand side of expression must be identifier"))
      }
    case a: Ast if a.isArrayOrMapLookup =>
      (evaluate(a.getAttribute("lhs")), evaluate(a.getAttribute("rhs"))) match {
        case (Success(a: WomArrayType), Success(WomIntegerType)) => Success(a.memberType)
        case (Success(m: WomMapType), Success(_: WomType)) => Success(m.valueType)
        case (Success(otherLhs), Success(_)) => Failure(new WomExpressionException(s"Invalid indexing target. You cannot index a value of type '${otherLhs.stableName}'"))
        case (f: Failure[_], _) => f
        case (_, f: Failure[_]) => f
      }
    case a: Ast if a.isFunctionCall =>
      val name = a.getAttribute("name").sourceString
      val params = a.params map evaluate
      functions.getFunction(name)(params)
  }

  def tupleAstToWdlType(a: Ast): Try[WomType] = {
    val unevaluatedElements = a.getAttribute("values").astListAsVector
    // Tuple 1 is equivalent to the value inside it. Enables nesting parens, e.g. (1 + 2) + 3
    if (unevaluatedElements.size == 1) {
      evaluate(unevaluatedElements.head)
    } else if (unevaluatedElements.size == 2) {
      for {
        left <- evaluate(unevaluatedElements.head)
        right <- evaluate(unevaluatedElements(1))
      } yield WomPairType(left, right)
    } else {
      Failure(new WomExpressionException(s"WDL does not currently support tuples with n > 2: ${a.toPrettyString}"))
    }
  }
}
