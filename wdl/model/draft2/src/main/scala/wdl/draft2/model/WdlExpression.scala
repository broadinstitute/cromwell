package wdl.draft2.model

import cats.data.Validated.Valid
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.draft2.model.AstTools.{EnhancedAstNode, VariableReference}
import wdl.draft2.model.WdlExpression._
import wdl.draft2.model.formatter.SyntaxHighlighter
import wdl.draft2.model.types.WdlExpressionType
import wdl.draft2.model.expression._
import wdl.draft2.model.formatter.NullSyntaxHighlighter
import wdl.draft2.model.{FullyQualifiedName => _}
import wdl.draft2.parser.WdlParser
import wdl.draft2.parser.WdlParser.{Ast, AstList, AstNode, Terminal}
import wom.core._
import wom.expression._
import wom.graph._
import wom.graph.expression.AnonymousExpressionNode
import wom.graph.expression.AnonymousExpressionNode.AnonymousExpressionConstructor
import wom.types.{WomAnyType, WomType}
import wom.values.{WomFile, WomFloat, WomValue}

import scala.collection.JavaConverters._
import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.language.postfixOps
import scala.util.{Failure, Try}

case object NoLookup extends ScopedLookupFunction {
  def apply(value: String): WomValue =
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
    def params: Vector[AstNode] = Option(ast.getAttribute("params")).map(_.asInstanceOf[AstList].asScala.toVector).getOrElse(Vector.empty)
    def name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    def isFunctionCallWithFirstParameterBeingFile = ast.isFunctionCall && ast.params.nonEmpty && WdlFunctionsWithFirstParameterBeingFile.contains(ast.functionName)
    def isGlobFunctionCall = ast.isFunctionCall && ast.params.size == 1 && "glob".equals(ast.functionName)
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

    def isTerminal: Boolean = astNode.isInstanceOf[Terminal]
  }

  val parser = new WdlParser()

  /** Maps from a locally qualified name to a WomValue. */
  type ScopedLookupFunction = String => WomValue

  val BinaryOperators = Set(
    "Add", "Subtract", "Multiply", "Divide", "Remainder",
    "GreaterThan", "LessThan", "GreaterThanOrEqual", "LessThanOrEqual",
    "Equals", "NotEquals", "LogicalAnd", "LogicalOr"
  )

  val UnaryOperators = Set("LogicalNot", "UnaryPlus", "UnaryNegation")

  val WdlFunctionsWithFirstParameterBeingFile: Seq[String] = Seq(
    "read_int",
    "read_string",
    "read_float",
    "read_boolean",
    "read_lines",
    "read_map",
    "read_object",
    "read_tsv",
    "read_json",
    "size"
  )

  def evaluate(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions[WomValue]): Try[WomValue] =
    ValueEvaluator(lookup, functions).evaluate(ast)

  def evaluateFiles(ast: AstNode, lookup: ScopedLookupFunction, functions: WdlFunctions[WomValue], coerceTo: WomType = WomAnyType) =
    FileEvaluator(expression.ValueEvaluator(lookup, functions), coerceTo).evaluate(ast)

  def evaluateType(ast: AstNode, lookup: (String) => WomType, functions: WdlFunctions[WomType], from: Option[Scope] = None) =
    TypeEvaluator(lookup, functions, from).evaluate(ast)

  def fromString(expression: WorkflowSource): WdlExpression = {
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
      case TernaryIf(condition, ifTrue, ifFalse) =>
        val c = toString(condition, highlighter)
        val t = toString(ifTrue, highlighter)
        val f = toString(ifFalse, highlighter)
        s"if $c then $t else $f"
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

case class WdlExpression(ast: AstNode) extends WomValue {
  override val womType: WomType = WdlExpressionType

  def evaluate(lookup: ScopedLookupFunction, functions: WdlFunctions[WomValue]): Try[WomValue] =
    WdlExpression.evaluate(ast, lookup, functions)

  def evaluateFiles(lookup: ScopedLookupFunction, functions: WdlFunctions[WomValue], coerceTo: WomType): Try[Seq[WomFile]] =
    WdlExpression.evaluateFiles(ast, lookup, functions, coerceTo)

  def evaluateType(lookup: (String) => WomType, functions: WdlFunctions[WomType], from: Option[Scope] = None): Try[WomType] =
    WdlExpression.evaluateType(ast, lookup, functions, from)

  def containsFunctionCall = ast.containsFunctionCalls

  def toString(highlighter: SyntaxHighlighter): String = {
    WdlExpression.toString(ast, highlighter)
  }

  override def toWomString: String = toString(NullSyntaxHighlighter)

  def prerequisiteCallNames: Set[FullyQualifiedName] = {
    this.topLevelMemberAccesses map { _.lhsString }
  }
  def topLevelMemberAccesses: Set[MemberAccess] = AstTools.findTopLevelMemberAccesses(ast) map { MemberAccess(_) } toSet
  def variableReferences(from: Scope): Iterable[VariableReference] = AstTools.findVariableReferences(ast, from)
}

/**
  *
  * @param wdlExpression The wrapped WdlExpression.
  * @param from The Scope in which the WdlExpression is found, needed to adjust member access expressions located in
  *             conditionals (wrapped in optionals) or scatters (wrapped in arrays).
  */
final case class WdlWomExpression(wdlExpression: WdlExpression, from: Scope) extends WomExpression {
  override def sourceString = wdlExpression.valueString
  override def inputs: Set[String] = wdlExpression.variableReferences(from) map { _.referencedVariableName } toSet

  override def evaluateValue(variableValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    lazy val wdlFunctions = WdlStandardLibraryFunctions.fromIoFunctionSet(ioFunctionSet)
    wdlExpression.evaluate(variableValues.apply, wdlFunctions).toErrorOr
  }

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] =
    // All current usages of WdlExpression#evaluateType trace back to WdlNamespace, but this is not the
    // case in the brave new WOM-world.
    wdlExpression.evaluateType(inputTypes.apply, new WdlStandardLibraryFunctionsType, Option(from)).toErrorOr

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {
    lazy val wdlFunctions = new WdlStandardLibraryFunctions {
      override def readFile(path: String): String = Await.result(ioFunctionSet.readFile(path, None, failOnOverflow = false), Duration.Inf)

      override def writeFile(path: String, content: String): Try[WomFile] = Try(Await.result(ioFunctionSet.writeFile(path, content), Duration.Inf))

      override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = ioFunctionSet.stdout(params)

      override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = ioFunctionSet.stderr(params)

      override def globHelper(pattern: String): Seq[String] = Await.result(ioFunctionSet.glob(pattern), Duration.Inf)

      override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = Failure(new Exception("You shouldn't call 'size' from a FileEvaluator"))
    }
    wdlExpression.evaluateFiles(inputTypes.apply, wdlFunctions, coerceTo).toErrorOr.map(_.toSet[WomFile])
  }
}

object WdlWomExpression {

  /**
    * Links up inputs to an expression ready for it to be instantiated.
    *
    * If the input is found in an outer scope, we also make a new input node in the inner graph to represent it.
    */
  def findInputsforExpression(expression: WdlWomExpression,
                              innerLookup: Map[String, GraphNodePort.OutputPort],
                              outerLookup: Map[String, GraphNodePort.OutputPort],
                              preserveIndexForOuterLookups: Boolean,
                              owningScope: Scope): ErrorOr[Map[String, GraphNodePort.OutputPort]] = {

    def resolveVariable(v: AstTools.VariableReference): ErrorOr[(String, GraphNodePort.OutputPort)] = {
      val name = v.referencedVariableName
      (innerLookup.get(name), outerLookup.get(name)) match {
        case (Some(port), _) =>
          // If we can find the value locally, use it.
          // It might be a local value or it might be an OGIN already created by another Node for an outerLookup.
          Valid(name -> port)
        case (None, Some(port)) => Valid(name -> OuterGraphInputNode(WomIdentifier(name), port, preserveIndexForOuterLookups).singleOutputPort)
        case (None, None) => s"No input $name found evaluating inputs for expression ${expression.wdlExpression.toWomString}".invalidNel
      }
    }

    for {
      resolvedVariables <- expression.wdlExpression.variableReferences(owningScope).toList traverse resolveVariable
    } yield resolvedVariables.toMap
  }

  def toAnonymousExpressionNode[T <: AnonymousExpressionNode](nodeIdentifier: WomIdentifier,
                                                              expression: WdlWomExpression,
                                                              innerLookup: Map[String, GraphNodePort.OutputPort],
                                                              outerLookup: Map[String, GraphNodePort.OutputPort],
                                                              preserveIndexForOuterLookups: Boolean,
                                                              owningScope: Scope,
                                                              constructor: AnonymousExpressionConstructor[T]): ErrorOr[T] = {
    import common.validation.ErrorOr.ShortCircuitingFlatMap

    findInputsforExpression(expression, innerLookup, outerLookup, preserveIndexForOuterLookups, owningScope) flatMap { resolvedVariables =>
      AnonymousExpressionNode.fromInputMapping(nodeIdentifier, expression, resolvedVariables, constructor)
    }
  }
}

object TernaryIf {
  def unapply(arg: Ast): Option[(AstNode, AstNode, AstNode)] = {
    if (arg.getName.equals("TernaryIf")) {
      Option((arg.getAttribute("cond"), arg.getAttribute("iftrue"), arg.getAttribute("iffalse")))
    } else {
      None
    }
  }
}
