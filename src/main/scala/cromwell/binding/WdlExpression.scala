package cromwell.binding

import cromwell.binding.formatter.{NullSyntaxHighlighter, SyntaxHighlighter}
import cromwell.binding.types.{WdlType, WdlArrayType, WdlExpressionType}
import cromwell.binding.values._
import cromwell.parser.WdlParser
import cromwell.parser.WdlParser.{Ast, AstList, AstNode, Terminal}

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

class WdlExpressionException(message: String = null, cause: Throwable = null) extends RuntimeException(message, cause)

object WdlExpression {

  implicit class AstForExpressions(val ast: Ast) extends AnyVal {
    def isFunctionCall: Boolean = ast.getName == "FunctionCall"
    def isBinaryOperator: Boolean = BinaryOperators.contains(ast.getName)
    def isUnaryOperator: Boolean = UnaryOperators.contains(ast.getName)
    def functionName: String = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    def isMemberAccess: Boolean = ast.getName == "MemberAccess"
    def isArrayLiteral: Boolean = ast.getName == "ArrayLiteral"
    def isAllowedInnerFunctionCall: Boolean = ast.isFunctionCall && AllowedInnerFunctionCalls.contains(ast.functionName)
    def params = ast.getAttribute("params").asInstanceOf[AstList].asScala.toVector
    def name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString

    def isFunctionCallWithOneFileParameter: Boolean = (
      ast.isFunctionCall
        && ast.params.size == 1
        && PreEvaluableFunctionsWithSingleParameter.contains(ast.functionName))
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

  val PreEvaluableFunctionsWithSingleParameter: Seq[String] = Seq("read_int", "read_string")

  // Function calls which can be called from within other function calls.
  // e.g. read_int(stdout()) is fine but read_int(read_int(stdout)) is not.
  val AllowedInnerFunctionCalls = Seq[String]("stdout", "stderr")

  /**
   * Look within the expression for filenames which aren't explicitly listed as outputs. 
   */
  def preevaluateExpressionForFilenames(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions): Try[Seq[WdlFile]] = {
    ast match {
      // This is the only case which actually pre-evaluates anything. The other cases are just buck-passers:
      case a: Ast if a.isFunctionCallWithOneFileParameter =>
        val innerExpression = a.params.head
        innerExpression match {
          case innerAst: Ast if innerAst.containsFunctionCalls =>
            // Only allow function calls if there's only one and its on the allowed list:
            if (innerAst.isAllowedInnerFunctionCall) {
              Success(Seq.empty)
            } else {
              Failure(new IllegalArgumentException(s"Invalid inner function call in $innerExpression"))
            }
          case _ => evaluate(innerExpression, lookup, functions) match {
            case Success(value) => Success(Seq(WdlFile(value.valueString)))
            case Failure(e) => Failure(e)
          }
        }
      // Binary operators - find filenames in sub-expressions and merge the lists:
      case a: Ast if a.isBinaryOperator =>
        val lhs = preevaluateExpressionForFilenames(a.getAttribute("lhs"), lookup, functions)
        val rhs = preevaluateExpressionForFilenames(a.getAttribute("rhs"), lookup, functions)
        // Recurse both sides and add the lists together:
        lhs match {
          case Success(lhsValue) =>
            rhs match {
              case Success(rhsValue) => Success(Seq(lhsValue, rhsValue).flatten)
              case Failure(error) => Failure(error)
            }
          case Failure(error) => Failure(error)
        }
      case a: Ast if a.isUnaryOperator => preevaluateExpressionForFilenames(a.getAttribute("expression"), lookup, functions)

      case _ => Success(Seq())
    }
  }

  private def replaceInterpolationTag(string: String, tag: String, lookup: ScopedLookupFunction) =
    string.replace(tag, lookup(tag.substring(2, tag.length - 1)).valueString)

  def interpolate(str: String, lookup: ScopedLookupFunction): String =
    "\\$\\{([a-zA-Z]([a-zA-Z0-9_])*)\\}".r.findAllIn(str).foldLeft(str) {replaceInterpolationTag(_, _, lookup)}

  def evaluate(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions, interpolateStrings: Boolean = false): Try[WdlValue] = {
    ast match {
      case t: Terminal if t.getTerminalStr == "identifier" => Success(lookup(t.getSourceString))
      case t: Terminal if t.getTerminalStr == "integer" => Success(WdlInteger(t.getSourceString.toInt))
      case t: Terminal if t.getTerminalStr == "float" => Success(WdlFloat(t.getSourceString.toDouble))
      case t: Terminal if t.getTerminalStr == "boolean" => Success(WdlBoolean(t.getSourceString == "true"))
      case t: Terminal if t.getTerminalStr == "string" =>
        val strValue = if (interpolateStrings) interpolate(t.getSourceString, lookup) else t.getSourceString
        Success(WdlString(strValue))
      case a: Ast if a.isBinaryOperator =>
        val lhs = evaluate(a.getAttribute("lhs"), lookup, functions, interpolateStrings)
        val rhs = evaluate(a.getAttribute("rhs"), lookup, functions, interpolateStrings)
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
        val expression = evaluate(a.getAttribute("expression"), lookup, functions, interpolateStrings)
        a.getName match {
          case "LogicalNot" => for(e <- expression) yield e.not.get
          case "UnaryPlus" => for(e <- expression) yield e.unaryPlus.get
          case "UnaryNegation" => for(e <- expression) yield e.unaryMinus.get
          case _ => Failure(new WdlExpressionException(s"Invalid operator: ${a.getName}"))
        }
      case a: Ast if a.isArrayLiteral =>
        val evaluatedElements = a.getAttribute("values").asInstanceOf[AstList].asScala.toVector map {x =>
          evaluate(x, lookup, functions, interpolateStrings)
        }
        evaluatedElements.partition {_.isSuccess} match {
          case (_, failures) if failures.nonEmpty =>
            val message = failures.collect {case f: Failure[_] => f.exception.getMessage}.mkString("\n")
            Failure(new WdlExpressionException(s"Could not evaluate expression:\n$message"))
          case (successes, _) =>
            successes.map{_.get.wdlType}.toSet match {
              case s:Set[WdlType] if s.isEmpty =>
                Failure(new WdlExpressionException(s"Can't have empty array declarations (can't infer type)"))
              case s:Set[WdlType] if s.size == 1 =>
                Success(WdlArray(WdlArrayType(s.head), successes.map{_.get}.toSeq))
              case _ =>
                Failure(new WdlExpressionException("Arrays must have homogeneous types"))
            }
        }
      case a: Ast if a.isMemberAccess =>
        a.getAttribute("rhs") match {
          case rhs:Terminal if rhs.getTerminalStr == "identifier" =>
            evaluate(a.getAttribute("lhs"), lookup, functions, interpolateStrings).flatMap {
              case o: WdlObject =>
                o.value.get(rhs.getSourceString) match {
                  case Some(v:WdlValue) => Success(v)
                  case None => Failure(new WdlExpressionException(s"Could not find key ${rhs.getSourceString}"))
                }
              case ns: WdlNamespace => Success(lookup(ns.importedAs.map {n => s"$n.${rhs.getSourceString}"}.getOrElse(rhs.getSourceString)))
              case _ => Failure(new WdlExpressionException("Left-hand side of expression must be a WdlObject or Namespace"))
            }
          case _ => Failure(new WdlExpressionException("Right-hand side of expression must be identifier"))
        }
      case a: Ast if a.isFunctionCall =>
        val name = a.getAttribute("name").asInstanceOf[Terminal].getSourceString
        val params = a.params map { evaluate(_, lookup, functions, interpolateStrings) }
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
      case a: Ast if a.isFunctionCall =>
        val params = a.params map { a => toString(a, highlighter) }
        s"${highlighter.function(a.name)}(${params.mkString(", ")})"
      case a: Ast if a.isMemberAccess =>
        val lhs = toString(a.getAttribute("lhs"), highlighter)
        val rhs = toString(a.getAttribute("rhs"), highlighter)
        s"$lhs.$rhs"
    }
  }
}

case class WdlExpression(ast: AstNode) extends WdlValue {

  import WdlExpression._

  override val wdlType = WdlExpressionType
  def evaluate(lookup: ScopedLookupFunction, functions: WdlFunctions, interpolateStrings: Boolean = false): Try[WdlValue] =
    WdlExpression.evaluate(ast, lookup, functions, interpolateStrings)
  def preevaluateExpressionForFilenames(lookup: ScopedLookupFunction, functions: WdlFunctions): Try[Seq[WdlFile]] =
    WdlExpression.preevaluateExpressionForFilenames(ast, lookup: ScopedLookupFunction, functions: WdlFunctions)
  def containsFunctionCall = ast.containsFunctionCalls
  def toString(highlighter: SyntaxHighlighter): String = {
    WdlExpression.toString(ast, highlighter)
  }
  override def toWdlString: String = toString(NullSyntaxHighlighter)
}

trait WdlFunctions {
  type WdlFunction = Seq[Try[WdlValue]] => Try[WdlValue]

  def getFunction(name: String): WdlFunction
}