package wdl.transforms.base.ast2wdlom

import cats.data.NonEmptyList
import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr.{ErrorOr, _}
import common.validation.Validation._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements._
import wom.values._

import scala.util.Try

object AstNodeToExpressionElement {

  type EngineFunctionMaker = Vector[ExpressionElement] => ErrorOr[ExpressionElement]

  def astNodeToExpressionElement(customEngineFunctionMakers: Map[String, EngineFunctionMaker]): CheckedAtoB[GenericAstNode, ExpressionElement] = {
    CheckedAtoB.fromErrorOr("parse expression")(convert(customEngineFunctionMakers) _)
  }

  protected def convert(customEngineFunctionMakers: Map[String, EngineFunctionMaker])(ast: GenericAstNode): ErrorOr[ExpressionElement] = {
    implicit val recursiveConverter = CheckedAtoB.fromErrorOr(convert(customEngineFunctionMakers) _)
    implicit val recursiveKvConverter = AstNodeToKvPair.astNodeToKvPair
    ast match {

      case t: GenericTerminal if asPrimitive.isDefinedAt((t.getTerminalStr, t.getSourceString)) => asPrimitive((t.getTerminalStr, t.getSourceString)).map(PrimitiveLiteralExpressionElement)
      case t: GenericTerminal if t.getTerminalStr == "identifier" => IdentifierLookup(t.getSourceString).validNel
      case t: GenericTerminal if t.getTerminalStr == "none" => NoneLiteralElement.validNel

      case a: GenericAst if a.getName == "StringLiteral" => handleStringLiteral(a)
      case a: GenericAst if lhsRhsOperators.contains(a.getName) => useValidatedLhsAndRhs(a, lhsRhsOperators(a.getName))
      case a: GenericAst if unaryOperators.contains(a.getName) => a.getAttributeAs[ExpressionElement]("expression").map(unaryOperators(a.getName)).toValidated
      case a: GenericAst if a.getName == "TupleLiteral" => (a.getAttributeAsVector[ExpressionElement]("values") flatMap {
        case pair if pair.length == 2 => PairLiteral(pair.head, pair(1)).validNelCheck
        case singleton if singleton.length == 1 => singleton.head.validNelCheck
        case more => s"No WDL support for ${more.size}-tuples in draft 3".invalidNelCheck
      }).toValidated
      case a: GenericAst if a.getName == "ArrayLiteral" => a.getAttributeAsVector[ExpressionElement]("values").toValidated.map(ArrayLiteral)
      case a: GenericAst if a.getName == "ArrayOrMapLookup" => {
        (a.getAttributeAs[ExpressionElement]("lhs").toValidated: ErrorOr[ExpressionElement],
          a.getAttributeAs[ExpressionElement]("rhs").toValidated: ErrorOr[ExpressionElement]) mapN IndexAccess
      }
      case a: GenericAst if a.getName == "MemberAccess" => handleMemberAccess(a)
      case a: GenericAst if a.getName == "ObjectLiteral" =>
        (for {
          objectKvs <- a.getAttributeAsVector[KvPair]("map")
          asMap = objectKvs.map(kv => kv.key -> kv.value).toMap
        } yield ObjectLiteral(asMap)).toValidated
      case a: GenericAst if a.getName == "MapLiteral" =>
        final case class MapKvPair(key: ExpressionElement, value: ExpressionElement)
        def convertOnePair(astNode: GenericAstNode): ErrorOr[MapKvPair] = astNode match {
          case a: GenericAst if a.getName == "ObjectKV" || a.getName == "MapLiteralKv" =>
            val keyValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("key").toValidated
            val valueValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("value").toValidated

            (keyValidation, valueValidation) mapN { (key, value) => MapKvPair(key, value) }
        }

        val astNodeToKvPair: CheckedAtoB[GenericAstNode, MapKvPair] = CheckedAtoB.fromErrorOr(convertOnePair)
        (for {
          mapKvs <- a.getAttributeAsVector[MapKvPair]("map")(astNodeToKvPair)
          asMap = mapKvs.map(kv => kv.key -> kv.value).toMap
        } yield MapLiteral(asMap)).toValidated
      case a: GenericAst if a.getName == "TernaryIf" =>
        val conditionValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("cond").toValidated
        val ifTrueValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("iftrue").toValidated
        val ifFalseValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("iffalse").toValidated
        (conditionValidation, ifTrueValidation, ifFalseValidation) mapN { (cond, ifTrue, ifFalse) => TernaryIf(cond, ifTrue, ifFalse) }
      case a: GenericAst if a.getName == "FunctionCall" =>
        val functionNameValidation: ErrorOr[String] = a.getAttributeAs[String]("name").toValidated
        val argsValidation: ErrorOr[Vector[ExpressionElement]] = a.getAttributeAsVector[ExpressionElement]("params").toValidated

        val allEngineFunctionMakers = engineFunctionMakers ++ customEngineFunctionMakers
        (functionNameValidation, argsValidation) flatMapN {
          case (name, params) if allEngineFunctionMakers.contains(name) => allEngineFunctionMakers(name).apply(params)
          case (other, _) => s"Unknown engine function: '$other'".invalidNel
        }


      case unknownTerminal: GenericTerminal => s"No rule available to create ExpressionElement from terminal: ${unknownTerminal.getTerminalStr} ${unknownTerminal.getSourceString}".invalidNel
      case unknownAst: GenericAst => s"No rule available to create ExpressionElement from Ast: ${unknownAst.getName}".invalidNel
    }
  }

  private type BinaryOperatorElementMaker = (ExpressionElement, ExpressionElement) => ExpressionElement

  private val lhsRhsOperators: Map[String, BinaryOperatorElementMaker] = Map(
    "LogicalOr" -> LogicalOr.apply,
    "LogicalAnd" -> LogicalAnd.apply,
    "Equals" -> Equals.apply,
    "NotEquals" -> NotEquals.apply,
    "LessThan" -> LessThan.apply,
    "LessThanOrEqual" -> LessThanOrEquals.apply,
    "GreaterThan" -> GreaterThan.apply,
    "GreaterThanOrEqual" -> GreaterThanOrEquals.apply,
    "Add" -> Add.apply,
    "Subtract" -> Subtract.apply,
    "Multiply" -> Multiply.apply,
    "Divide" -> Divide.apply,
    "Remainder" -> Remainder.apply
  )

  private def useValidatedLhsAndRhs(a: GenericAst, combiner: BinaryOperatorElementMaker)
                                   (implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]): ErrorOr[ExpressionElement] = {
    val lhsValidation: ErrorOr[ExpressionElement] =
      a.getAttributeAs[ExpressionElement]("lhs").toValidated.contextualizeErrors(s"read left hand side of ${a.getName} expression")
    val rhsValidation: ErrorOr[ExpressionElement] =
      a.getAttributeAs[ExpressionElement]("rhs").toValidated.contextualizeErrors(s"read right-hand side of ${a.getName} expression")

    (lhsValidation, rhsValidation) mapN { combiner }
  }

  private type UnaryOperatorElementMaker = ExpressionElement => ExpressionElement

  private val unaryOperators: Map[String, UnaryOperatorElementMaker] = Map(
    "LogicalNot" -> LogicalNot.apply,
    "UnaryPlus" -> UnaryPlus.apply,
    "UnaryNegation" -> UnaryNegation.apply
  )

  private val asPrimitive: PartialFunction[(String, String), ErrorOr[WomPrimitive]] = {
    case ("integer", i) => Try(WomInteger(i.toInt)).toErrorOr
    case ("float", f) => Try(WomFloat(f.toDouble)).toErrorOr
    case ("boolean", b) => Try(WomBoolean(b.toBoolean)).toErrorOr
    case ("string", s) => WomString(s).validNel
  }

  private val engineFunctionMakers: Map[String, Vector[ExpressionElement] => ErrorOr[ExpressionElement]] = Map(
    // 0-param functions:
    "stdout" -> validateNoParamEngineFunction(StdoutElement, "stdout"),
    "stderr" -> validateNoParamEngineFunction(StderrElement, "stderr"),

    // 1-param functions:
    "read_lines" -> validateOneParamEngineFunction(ReadLines, "read_lines"),
    "read_tsv" -> validateOneParamEngineFunction(ReadTsv, "read_tsv"),
    "read_map" -> validateOneParamEngineFunction(ReadMap, "read_map"),
    "read_object" -> validateOneParamEngineFunction(ReadObject, "read_object"),
    "read_objects" -> validateOneParamEngineFunction(ReadObjects, "read_objects"),
    "read_json" -> validateOneParamEngineFunction(ReadJson, "read_json"),
    "read_int" -> validateOneParamEngineFunction(ReadInt, "read_int"),
    "read_string" -> validateOneParamEngineFunction(ReadString, "read_string"),
    "read_float" -> validateOneParamEngineFunction(ReadFloat, "read_float"),
    "read_boolean" -> validateOneParamEngineFunction(ReadBoolean, "read_boolean"),
    "write_lines" -> validateOneParamEngineFunction(WriteLines, "write_lines"),
    "write_tsv" -> validateOneParamEngineFunction(WriteTsv, "write_tsv"),
    "write_map" -> validateOneParamEngineFunction(WriteMap, "write_map"),
    "write_object" -> validateOneParamEngineFunction(WriteObject, "write_object"),
    "write_objects" -> validateOneParamEngineFunction(WriteObjects, "write_objects"),
    "write_json" -> validateOneParamEngineFunction(WriteJson, "write_json"),
    "range" -> validateOneParamEngineFunction(Range, "range"),
    "transpose" -> validateOneParamEngineFunction(Transpose, "transpose"),
    "length" -> validateOneParamEngineFunction(Length, "length"),
    "flatten" -> validateOneParamEngineFunction(Flatten, "flatten"),
    "select_first" -> validateOneParamEngineFunction(SelectFirst, "select_first"),
    "select_all" -> validateOneParamEngineFunction(SelectAll, "select_all"),
    "defined" -> validateOneParamEngineFunction(Defined, "defined"),
    "floor" -> validateOneParamEngineFunction(Floor, "floor"),
    "ceil" -> validateOneParamEngineFunction(Ceil, "ceil"),
    "round" -> validateOneParamEngineFunction(Round, "round"),
    "glob" -> validateOneParamEngineFunction(Glob, "glob"),

    // 1- or 2-param functions:
    "size" -> validateOneOrTwoParamEngineFunction(Size, "size"),
    "basename" -> validateOneOrTwoParamEngineFunction(Basename, "basename"),

    // 2-param functions:
    "zip" -> validateTwoParamEngineFunction(Zip, "zip"),
    "cross" -> validateTwoParamEngineFunction(Cross, "cross"),
    "prefix" -> validateTwoParamEngineFunction(Prefix, "prefix"),

    // 3-param functions:
    "sub" -> validateThreeParamEngineFunction(Sub, "sub")
  )

  private def validateNoParamEngineFunction(element: ExpressionElement, functionName: String)
                                           (params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] =
    if (params.isEmpty) {
      element.validNel
    } else {
      s"Function $functionName expects 0 arguments but got ${params.size}".invalidNel
    }

  def validateOneParamEngineFunction(elementMaker: ExpressionElement => ExpressionElement, functionName: String)
                                            (params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] =
    if (params.size == 1) {
      elementMaker.apply(params.head).validNel
    } else {
      s"Function $functionName expects exactly 1 argument but got ${params.size}".invalidNel
    }

  private def validateOneOrTwoParamEngineFunction(elementMaker: (ExpressionElement, Option[ExpressionElement]) => ExpressionElement,
                                                  functionName: String)
                                                 (params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] =
    if (params.size == 1) {
      elementMaker(params.head, None).validNel
    } else if (params.size == 2) {
      elementMaker(params.head, Option(params(1))).validNel
    } else {
      s"Function $functionName expects 1 or 2 arguments but got ${params.size}".invalidNel
    }

  private def validateTwoParamEngineFunction(elementMaker: (ExpressionElement, ExpressionElement) => ExpressionElement, functionName: String)
                                            (params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] =
    if (params.size == 2) {
      elementMaker.apply(params.head, params(1)).validNel
    } else {
      s"Function $functionName expects exactly 2 arguments but got ${params.size}".invalidNel
    }

  private def validateThreeParamEngineFunction(elementMaker: (ExpressionElement, ExpressionElement, ExpressionElement) => ExpressionElement, functionName: String)
                                              (params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] =
    if (params.size == 3) {
      elementMaker.apply(params.head, params(1), params(2)).validNel
    } else {
      s"Function $functionName expects exactly 3 arguments but got ${params.size}".invalidNel
    }

  private def handleMemberAccess(ast: GenericAst)
                                (implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]): ErrorOr[ExpressionElement] = {

    // Internal simplify method:
    // If the left-hand side is another member access, we can simplify them together.
    // If not, we can build a new member access:
    def simplify(leftExpression: ExpressionElement, suffix: String): ExpressionElement = leftExpression match {
      case ExpressionMemberAccess(lefterExpression, tail) => ExpressionMemberAccess(lefterExpression, NonEmptyList(tail.head, tail.tail :+ suffix))
      case IdentifierMemberAccess(first, second, tail) => IdentifierMemberAccess(first, second, tail :+ suffix)
      case IdentifierLookup(identifier) => IdentifierMemberAccess(identifier, suffix, Vector.empty)
      case _ => ExpressionMemberAccess(leftExpression, NonEmptyList(suffix, List.empty))
    }

    val leftValidation: ErrorOr[ExpressionElement] = ast.getAttributeAs[ExpressionElement]("value").toValidated
    val rightValidation: ErrorOr[String] = ast.getAttributeAs[String]("member").toValidated

    (leftValidation, rightValidation) mapN { (lhs, rhs) =>
      simplify(lhs, rhs)
    }

  }

  private def handleStringLiteral(ast: GenericAst)
                                 (implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]): ErrorOr[ExpressionElement] = {

    implicit val astNodeToStringPiece: CheckedAtoB[GenericAstNode, StringPiece] = AstNodeToStringPiece.astNodeToStringPiece(Some(astNodeToExpressionElement))
    ast.getAttributeAsVector[StringPiece]("pieces").toValidated map { pieces =>
      if (pieces.isEmpty) {
        StringLiteral("")
      } else if (pieces.size == 1) {
        pieces.head match {
          case s: StringLiteral => s
          case e: StringEscapeSequence => StringLiteral(e.unescape)
          case _ => StringExpression(pieces)
        }
      } else {
        StringExpression(pieces)
      }
    }
  }
}
