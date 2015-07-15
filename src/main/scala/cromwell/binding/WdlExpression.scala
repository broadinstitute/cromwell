package cromwell.binding

import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.formatter.{NullSyntaxHighlighter, SyntaxHighlighter}
import cromwell.binding.types.WdlExpressionType
import cromwell.binding.values._
import cromwell.parser.WdlParser
import cromwell.parser.WdlParser.{Ast, AstList, AstNode, Terminal}

import scala.collection.JavaConverters._
import scala.util.{Success, Try}

class WdlExpressionException(message: String = null, cause: Throwable = null) extends RuntimeException(message, cause)

object WdlExpression {
  val parser = new WdlParser()

  /** Maps from a locally qualified name to a WdlValue. */
  type ScopedLookupFunction = String => WdlValue

  def binaryOperators = Set(
    "Add", "Subtract", "Multiply", "Divide", "Remainder",
    "GreaterThan", "LessThan", "GreaterThanOrEqual", "LessThanOrEqual",
    "Equals", "NotEquals", "LogicalAnd", "LogicalOr"
  )

  val unaryOperators = Set("LogicalNot", "UnaryPlus", "UnaryNegation")

  def evaluate(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions): Try[WdlValue] = {
    ast match {
      case t: Terminal if t.getTerminalStr == "identifier" => Success(lookup(t.getSourceString))
      case t: Terminal if t.getTerminalStr == "integer" => Success(WdlInteger(t.getSourceString.toInt))
      case t: Terminal if t.getTerminalStr == "float" => Success(WdlFloat(t.getSourceString.toDouble))
      case t: Terminal if t.getTerminalStr == "boolean" => Success(WdlBoolean(t.getSourceString == "true"))
      case t: Terminal if t.getTerminalStr == "string" => Success(WdlString(t.getSourceString))
      case a: Ast if binaryOperators.contains(a.getName) =>
        val lhs = evaluate(a.getAttribute("lhs"), lookup, functions)
        val rhs = evaluate(a.getAttribute("rhs"), lookup, functions)
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
          case _ => throw new WdlExpressionException(s"Invalid operator: ${a.getName}")
        }
      case a: Ast if unaryOperators.contains(a.getName) =>
        val expression = evaluate(a.getAttribute("expression"), lookup, functions)
        a.getName match {
          case "LogicalNot" => for(e <- expression) yield e.not.get
          case "UnaryPlus" => for(e <- expression) yield e.unaryPlus.get
          case "UnaryNegation" => for(e <- expression) yield e.unaryMinus.get
          case _ => throw new WdlExpressionException(s"Invalid operator: ${a.getName}")
        }
      case a: Ast if a.getName == "MemberAccess" =>
        val rhs = a.getAttribute("rhs") match {
          case x:Terminal if x.getTerminalStr == "identifier" => x.getSourceString
          case _ => throw new WdlExpressionException("Right-hand side of expression must be identifier")
        }
        evaluate(a.getAttribute("lhs"), lookup, functions).map {
          case o: WdlObject => o.value.getOrElse(rhs, throw new WdlExpressionException(s"Could not find key $rhs"))
          case ns: WdlNamespace => lookup(ns.namespace.map {n => s"$n.$rhs"}.getOrElse(rhs))
          case _ => throw new WdlExpressionException("Left-hand side of expression must be a WdlObject or Namespace")
        }
      case a: Ast if a.getName == "FunctionCall" =>
        val name = a.getAttribute("name").asInstanceOf[Terminal].getSourceString
        val params = a.getAttribute("params").asInstanceOf[AstList].asScala.toVector map {
          evaluate(_, lookup, functions)
        }
        functions.getFunction(name)(params)
    }
  }

  def fromString(expression: WdlSource): WdlExpression = {
    val tokens = parser.lex(expression, "string")
    val terminalMap = (tokens.asScala.toVector map {(_, expression)}).toMap
    val parseTree = parser.parse_e(tokens, new WdlSyntaxErrorFormatter(terminalMap))
    new WdlExpression(parseTree.toAst)
  }

  def toString(ast: AstNode, highlighter: SyntaxHighlighter = NullSyntaxHighlighter): String = {
    ast match {
      case t: Terminal if t.getTerminalStr == "identifier" => t.getSourceString
      case t: Terminal if t.getTerminalStr == "integer" => t.getSourceString
      case t: Terminal if t.getTerminalStr == "float" => t.getSourceString
      case t: Terminal if t.getTerminalStr == "string" => s""""${t.getSourceString}""""
      case a:Ast if binaryOperators.contains(a.getName) => {
        val lhs = toString(a.getAttribute("lhs"), highlighter)
        val rhs = toString(a.getAttribute("rhs"), highlighter)
        a.getName match {
          case "Add" => s"$lhs + $rhs"
          case "Subtract" => s"$lhs - $rhs"
          case "Multiply" => s"$lhs * $rhs"
          case "Divide" => s"$lhs / $rhs"
          case "Remainder" => s"$lhs % $rhs"
          case "Equals" => s"$lhs == $rhs"
          case "NotEquals" => s"$lhs != $rhs"
          case "LessThan" => s"$lhs < $rhs"
          case "LessThanOrEqual" => s"$lhs <= $rhs"
          case "GreaterThan" => s"$lhs > $rhs"
          case "GreaterThanOrEqual" => s"$lhs >= $rhs"
          case "LogicalOr" => s"$lhs || $rhs"
          case "LogicalAnd" => s"$lhs && $rhs"
        }
      }
      case a:Ast if a.getName == "FunctionCall" => {
        val name = a.getAttribute("name").asInstanceOf[Terminal].getSourceString
        val params = a.getAttribute("params").asInstanceOf[AstList].asScala.toVector.map {a => toString(a, highlighter)}
        s"${highlighter.function(name)}(${params.mkString(", ")})"
      }
      case a:Ast if a.getName == "MemberAccess" => {
        val lhs = toString(a.getAttribute("lhs"), highlighter)
        val rhs = toString(a.getAttribute("rhs"), highlighter)
        s"$lhs.$rhs"
      }
    }
  }
}

case class WdlExpression(ast: AstNode) extends WdlValue {
  override val wdlType = WdlExpressionType
  def evaluate(lookup: ScopedLookupFunction, functions: WdlFunctions): Try[WdlValue] =
    WdlExpression.evaluate(ast, lookup, functions)
  def toString(highlighter: SyntaxHighlighter): String = {
    WdlExpression.toString(ast, highlighter)
  }
  override def toWdlString: String = toString(NullSyntaxHighlighter)
}

trait WdlFunctions {
  type WdlFunction = Seq[Try[WdlValue]] => Try[WdlValue]

  def getFunction(name: String): WdlFunction
}
