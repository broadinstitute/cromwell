package wdl4s

import lenthall.exception.ThrowableAggregation
import wdl4s.AstTools.EnhancedAstNode
import wdl4s.WdlExpression._
import wdl4s.expression._
import wdl4s.formatter.{NullSyntaxHighlighter, SyntaxHighlighter}
import wdl4s.parser.WdlParser
import wdl4s.parser.WdlParser.{Ast, AstList, AstNode, Terminal}
import wdl4s.types._
import wdl4s.values._

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.Try

class WdlExpressionException(message: String = null, cause: Throwable = null) extends RuntimeException(message, cause)

case object NoLookup extends ScopedLookupFunction {
  def apply(value: String): WdlValue =
    throw new UnsupportedOperationException(s"No identifiers should be looked up: $value")
}

object WdlExpression {

  implicit class AstForExpressions(val ast: Ast) extends AnyVal {
    def isFunctionCall: Boolean = ast.getName == "FunctionCall"
    def isBinaryOperator: Boolean = BinaryOperators.contains(ast.getName)
    def isUnaryOperator: Boolean = UnaryOperators.contains(ast.getName)
    def functionName: String = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    def isMemberAccess: Boolean = ast.getName == "MemberAccess"
    def isArrayLiteral: Boolean = ast.getName == "ArrayLiteral"
    def isTupleLiteral: Boolean = ast.getName == "TupleLiteral"
    def isMapLiteral: Boolean = ast.getName == "MapLiteral"
    def isObjectLiteral: Boolean = ast.getName == "ObjectLiteral"
    def isArrayOrMapLookup: Boolean = ast.getName == "ArrayOrMapLookup"
    def params = ast.getAttribute("params").asInstanceOf[AstList].asScala.toVector
    def name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    def isFunctionCallWithOneParameter = ast.isFunctionCall && ast.params.size == 1
    def isFunctionCallWithOneFileParameter = isFunctionCallWithOneParameter && WdlFunctionsWithSingleFileParameter.contains(ast.functionName)
    def isGlobFunctionCall = isFunctionCallWithOneParameter && "glob".equals(ast.functionName)
  }

  implicit class AstNodeForExpressions(val astNode: AstNode) extends AnyVal {
    def containsFunctionCalls: Boolean =
      astNode match {
        case a: Ast if a.isFunctionCall => true
        case a: Ast if a.isBinaryOperator =>
          val lhs = a.getAttribute("lhs")
          val rhs = a.getAttribute("rhs")
          lhs.containsFunctionCalls || rhs.containsFunctionCalls
        case a: Ast if a.isUnaryOperator =>
          val rhs = a.getAttribute("expression")
          rhs.containsFunctionCalls
        case _ => false
      }
  }

  val parser = new WdlParser()

  /** Maps from a locally qualified name to a WdlValue. */
  type ScopedLookupFunction = String => WdlValue

  val BinaryOperators = Set(
    "Add", "Subtract", "Multiply", "Divide", "Remainder",
    "GreaterThan", "LessThan", "GreaterThanOrEqual", "LessThanOrEqual",
    "Equals", "NotEquals", "LogicalAnd", "LogicalOr"
  )

  val UnaryOperators = Set("LogicalNot", "UnaryPlus", "UnaryNegation")

  val WdlFunctionsWithSingleFileParameter: Seq[String] = Seq(
    "read_int",
    "read_string",
    "read_float",
    "read_boolean",
    "read_lines",
    "read_map",
    "read_object",
    "read_tsv",
    "size"
  )

  def evaluate(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions[WdlValue]): Try[WdlValue] =
    ValueEvaluator(lookup, functions).evaluate(ast)

  def evaluateFiles(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions[WdlValue], coerceTo: WdlType = WdlAnyType) =
    FileEvaluator(ValueEvaluator(lookup, functions), coerceTo).evaluate(ast)

  def evaluateType(ast: AstNode, lookup: (String) => WdlType, functions: WdlFunctions[WdlType], from: Option[Scope] = None) =
    TypeEvaluator(lookup, functions, from).evaluate(ast)

  def fromString(expression: WdlSource): WdlExpression = {
    val tokens = parser.lex(expression, "string")
    val terminalMap = (tokens.asScala.toVector map {(_, expression)}).toMap
    val parseTree = parser.parse_e(tokens, WdlSyntaxErrorFormatter(terminalMap))
    new WdlExpression(parseTree.toAst)
  }

  def toString(ast: AstNode, highlighter: SyntaxHighlighter = NullSyntaxHighlighter): String = {
    ast match {
      case t: Terminal if Seq("identifier", "integer", "float", "boolean").contains(t.getTerminalStr) => t.getSourceString
      case t: Terminal if t.getTerminalStr == "string" => s""""${t.getSourceString.replaceAll("\"", "\\" + "\"")}""""
      case a:Ast if a.isBinaryOperator =>
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
      case a: Ast if a.isUnaryOperator =>
        val expression = toString(a.getAttribute("expression"), highlighter)
        a.getName match {
          case "LogicalNot" => s"!$expression"
          case "UnaryPlus" => s"+$expression"
          case "UnaryNegation" => s"-$expression"
        }
      case a: Ast if a.isArrayLiteral =>
        val evaluatedElements = a.getAttribute("values").astListAsVector map {x => toString(x, highlighter)}
        s"[${evaluatedElements.mkString(",")}]"
      case a: Ast if a.isTupleLiteral =>
        val evaluatedElements = a.getAttribute("values").astListAsVector map { x => toString(x, highlighter)}
        s"(${evaluatedElements.mkString(", ")})"
      case a: Ast if a.isMapLiteral =>
        val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
          val key = toString(kv.asInstanceOf[Ast].getAttribute("key"), highlighter)
          val value = toString(kv.asInstanceOf[Ast].getAttribute("value"), highlighter)
          s"$key:$value"
        }
        s"{${evaluatedMap.mkString(",")}}"
      case a: Ast if a.isMemberAccess =>
        val lhs = toString(a.getAttribute("lhs"), highlighter)
        val rhs = toString(a.getAttribute("rhs"), highlighter)
        s"$lhs.$rhs"
      case a: Ast if a.isArrayOrMapLookup =>
        val lhs = toString(a.getAttribute("lhs"), highlighter)
        val rhs = toString(a.getAttribute("rhs"), highlighter)
        s"$lhs[$rhs]"
      case a: Ast if a.isFunctionCall =>
        val params = a.params map { a => toString(a, highlighter) }
        s"${highlighter.function(a.name)}(${params.mkString(", ")})"
    }
  }
}

case class WdlExpression(ast: AstNode) extends WdlValue {
  override val wdlType = WdlExpressionType

  def evaluate(lookup: ScopedLookupFunction, functions: WdlFunctions[WdlValue]): Try[WdlValue] =
    WdlExpression.evaluate(ast, lookup, functions)

  def evaluateFiles(lookup: ScopedLookupFunction, functions: WdlFunctions[WdlValue], coerceTo: WdlType): Try[Seq[WdlFile]] =
    WdlExpression.evaluateFiles(ast, lookup, functions, coerceTo)

  def evaluateType(lookup: (String) => WdlType, functions: WdlFunctions[WdlType], from: Option[Scope] = None): Try[WdlType] =
    WdlExpression.evaluateType(ast, lookup, functions, from)

  def containsFunctionCall = ast.containsFunctionCalls

  def toString(highlighter: SyntaxHighlighter): String = {
    WdlExpression.toString(ast, highlighter)
  }

  override def toWdlString: String = toString(NullSyntaxHighlighter)

  def prerequisiteCallNames: Set[FullyQualifiedName] = {
    this.topLevelMemberAccesses map { _.lhs }
  }
  def topLevelMemberAccesses: Set[MemberAccess] = AstTools.findTopLevelMemberAccesses(ast) map { MemberAccess(_) } toSet
  def variableReferences: Iterable[Terminal] = AstTools.findVariableReferences(ast)
}
