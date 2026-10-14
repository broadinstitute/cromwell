package wdl.transforms.v1_2.ast2wdlom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{
  AsMap,
  AsPairs,
  CollectByKey,
  Keys,
  Max,
  Min,
  Quote,
  ReadTsv,
  SelectFirst,
  Sep,
  SQuote,
  SubPosix,
  Suffix,
  Unzip,
  WriteTsv
}
import wdl.model.v1_2.elements.V1_2ExpressionElement._
import wdl.transforms.base.ast2wdlom.AstNodeToExpressionElement

object AstToNewExpressionElements {
  val newV1_2EngineFunctionMakers: Map[String, Vector[ExpressionElement] => ErrorOr[ExpressionElement]] = Map(
    // Functions new in WDL 1.1:
    "keys" -> AstNodeToExpressionElement.validateOneParamEngineFunction(Keys, "keys"),
    "as_map" -> AstNodeToExpressionElement.validateOneParamEngineFunction(AsMap, "as_map"),
    "as_pairs" -> AstNodeToExpressionElement.validateOneParamEngineFunction(AsPairs, "as_pairs"),
    "collect_by_key" -> AstNodeToExpressionElement.validateOneParamEngineFunction(CollectByKey, "collect_by_key"),
    "min" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Min, "min"),
    "max" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Max, "max"),
    "sep" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Sep, "sep"),
    "sub" -> AstNodeToExpressionElement.validateThreeParamEngineFunction(SubPosix, "sub"),
    "suffix" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Suffix, "suffix"),
    "quote" -> AstNodeToExpressionElement.validateOneParamEngineFunction(Quote, "quote"),
    "squote" -> AstNodeToExpressionElement.validateOneParamEngineFunction(SQuote, "squote"),
    "unzip" -> AstNodeToExpressionElement.validateOneParamEngineFunction(Unzip, "unzip"),

    // Functions new in WDL 1.2:
    "find" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Find, "find"),
    "matches" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Matches, "matches"),
    "join_paths" -> validateOneOrTwoParams(JoinPaths, "join_paths"),
    "contains" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Contains, "contains"),
    "chunk" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(Chunk, "chunk"),
    "contains_key" -> AstNodeToExpressionElement.validateTwoParamEngineFunction(ContainsKey, "contains_key"),
    "values" -> AstNodeToExpressionElement.validateOneParamEngineFunction(Values, "values"),

    // Functions with new overloads in WDL 1.2:
    "select_first" -> selectFirst,
    "read_tsv" -> readTsv,

    // Operators new in WDL 1.2, presented to the expression converter as function calls:
    V1_2GenericAstNode.ExponentiationFunctionName -> AstNodeToExpressionElement.validateTwoParamEngineFunction(
      Exponentiation,
      V1_2GenericAstNode.ExponentiationFunctionName
    ),
    "write_tsv" -> writeTsv
    // `read_object`, `read_objects`, `write_object` and `write_objects` are deprecated in WDL 1.2 together with the
    // `Object` type, but they are part of the language again (WDL 1.1 Cromwell removed them).
    // TODO: Suggestion: report a deprecation warning when they are used, once language factories can report warnings.
  )

  private def validateOneOrTwoParams(elementMaker: (ExpressionElement, Option[ExpressionElement]) => ExpressionElement,
                                     functionName: String
  )(params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] = params match {
    case Vector(first) => elementMaker(first, None).validNel
    case Vector(first, second) => elementMaker(first, Option(second)).validNel
    case _ => s"Function $functionName expects 1 or 2 arguments but got ${params.size}".invalidNel
  }

  // `X select_first(Array[X?]+)` or `X select_first(Array[X?], X)`
  private def selectFirst(params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] = params match {
    case Vector(array) => SelectFirst(array).validNel
    case Vector(array, default) => SelectFirstWithDefault(array, default).validNel
    case _ => s"Function select_first expects 1 or 2 arguments but got ${params.size}".invalidNel
  }

  // `File write_tsv(Array[Array[String]]|Array[Struct])` or `File write_tsv(Array, Boolean[, Array[String]])`
  private def writeTsv(params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] = params match {
    case Vector(array) => WriteTsv(array).validNel
    case Vector(array, header) => WriteTsvWithHeader(array, header, None).validNel
    case Vector(array, header, columnNames) => WriteTsvWithHeader(array, header, Option(columnNames)).validNel
    case _ => s"Function write_tsv expects 1, 2 or 3 arguments but got ${params.size}".invalidNel
  }

  // `Array[Array[String]] read_tsv(File)` or `Array[Object] read_tsv(File, Boolean[, Array[String]])`
  private def readTsv(params: Vector[ExpressionElement]): ErrorOr[ExpressionElement] = params match {
    case Vector(file) => ReadTsv(file).validNel
    case Vector(file, header) => ReadTsvWithHeader(file, header, None).validNel
    case Vector(file, header, fieldNames) => ReadTsvWithHeader(file, header, Option(fieldNames)).validNel
    case _ => s"Function read_tsv expects 1, 2 or 3 arguments but got ${params.size}".invalidNel
  }
}
