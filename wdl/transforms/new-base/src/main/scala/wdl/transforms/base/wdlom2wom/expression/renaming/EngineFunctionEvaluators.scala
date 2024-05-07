package wdl.transforms.base.wdlom2wom.expression.renaming

import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements.ExpressionElement

object EngineFunctionEvaluators {

  implicit val stdoutRenamer: IdentifierLookupRenamer[StdoutElement.type] = forZeroParamFunction[StdoutElement.type]
  implicit val stderrRenamer: IdentifierLookupRenamer[StderrElement.type] = forZeroParamFunction[StderrElement.type]

  implicit val readLinesRenamer: IdentifierLookupRenamer[ReadLines] = forOneParamFunction(ReadLines)
  implicit val readTsvRenamer: IdentifierLookupRenamer[ReadTsv] = forOneParamFunction(ReadTsv)
  implicit val readMapRenamer: IdentifierLookupRenamer[ReadMap] = forOneParamFunction(ReadMap)
  implicit val readObjectRenamer: IdentifierLookupRenamer[ReadObject] = forOneParamFunction(ReadObject)
  implicit val readObjectsRenamer: IdentifierLookupRenamer[ReadObjects] = forOneParamFunction(ReadObjects)
  implicit val readJsonRenamer: IdentifierLookupRenamer[ReadJson] = forOneParamFunction(ReadJson)
  implicit val readIntRenamer: IdentifierLookupRenamer[ReadInt] = forOneParamFunction(ReadInt)
  implicit val readStringRenamer: IdentifierLookupRenamer[ReadString] = forOneParamFunction(ReadString)
  implicit val readFloatRenamer: IdentifierLookupRenamer[ReadFloat] = forOneParamFunction(ReadFloat)
  implicit val readBooleanRenamer: IdentifierLookupRenamer[ReadBoolean] = forOneParamFunction(ReadBoolean)
  implicit val writeLinesRenamer: IdentifierLookupRenamer[WriteLines] = forOneParamFunction(WriteLines)
  implicit val writeTsvRenamer: IdentifierLookupRenamer[WriteTsv] = forOneParamFunction(WriteTsv)
  implicit val writeMapRenamer: IdentifierLookupRenamer[WriteMap] = forOneParamFunction(WriteMap)
  implicit val writeObjectRenamer: IdentifierLookupRenamer[WriteObject] = forOneParamFunction(WriteObject)
  implicit val writeObjectsRenamer: IdentifierLookupRenamer[WriteObjects] = forOneParamFunction(WriteObjects)
  implicit val writeJsonRenamer: IdentifierLookupRenamer[WriteJson] = forOneParamFunction(WriteJson)
  implicit val rangeRenamer: IdentifierLookupRenamer[Range] = forOneParamFunction(Range)
  implicit val transposeRenamer: IdentifierLookupRenamer[Transpose] = forOneParamFunction(Transpose)
  implicit val lengthRenamer: IdentifierLookupRenamer[Length] = forOneParamFunction(Length)
  implicit val flattenRenamer: IdentifierLookupRenamer[Flatten] = forOneParamFunction(Flatten)
  implicit val selectFirstRenamer: IdentifierLookupRenamer[SelectFirst] = forOneParamFunction(SelectFirst)
  implicit val selectAllRenamer: IdentifierLookupRenamer[SelectAll] = forOneParamFunction(SelectAll)
  implicit val definedRenamer: IdentifierLookupRenamer[Defined] = forOneParamFunction(Defined)
  implicit val floorRenamer: IdentifierLookupRenamer[Floor] = forOneParamFunction(Floor)
  implicit val ceilRenamer: IdentifierLookupRenamer[Ceil] = forOneParamFunction(Ceil)
  implicit val roundRenamer: IdentifierLookupRenamer[Round] = forOneParamFunction(Round)
  implicit val globRenamer: IdentifierLookupRenamer[Glob] = forOneParamFunction(Glob)

  implicit val basenameRenamer: IdentifierLookupRenamer[Basename] = forOneOrTwoParamFunction(Basename)
  implicit val sizeRenamer: IdentifierLookupRenamer[Size] = forOneOrTwoParamFunction(Size)

  implicit val zipRenamer: IdentifierLookupRenamer[Zip] = forTwoParamFunction(Zip)
  implicit val crossRenamer: IdentifierLookupRenamer[Cross] = forTwoParamFunction(Cross)
  implicit val prefixRenamer: IdentifierLookupRenamer[Prefix] = forTwoParamFunction(Prefix)

  implicit val subRenamer: IdentifierLookupRenamer[Sub] = forThreeParamFunction(Sub)

  private def forZeroParamFunction[A <: ExpressionElement]: IdentifierLookupRenamer[A] =
    new IdentifierLookupRenamer[A] {
      override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): A = a
    }

  private def forOneParamFunction[A <: OneParamFunctionCallElement](
    constructor: ExpressionElement => A
  ): IdentifierLookupRenamer[A] = new IdentifierLookupRenamer[A] {
    override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
    ): A =
      constructor.apply(
        expressionElementRenamer.renameIdentifiers(a.param, renamingMap)(expressionElementRenamer)
      )
  }

  private def forOneOrTwoParamFunction[A <: OneOrTwoParamFunctionCallElement](
    constructor: (ExpressionElement, Option[ExpressionElement]) => A
  ): IdentifierLookupRenamer[A] = new IdentifierLookupRenamer[A] {
    override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
    ): A =
      constructor.apply(
        expressionElementRenamer.renameIdentifiers(a.firstParam, renamingMap)(expressionElementRenamer),
        a.secondParam.map(expressionElementRenamer.renameIdentifiers(_, renamingMap)(expressionElementRenamer))
      )
  }

  private def forTwoParamFunction[A <: TwoParamFunctionCallElement](
    constructor: (ExpressionElement, ExpressionElement) => A
  ): IdentifierLookupRenamer[A] = new IdentifierLookupRenamer[A] {
    override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
    ): A =
      constructor.apply(
        expressionElementRenamer.renameIdentifiers(a.arg1, renamingMap)(expressionElementRenamer),
        expressionElementRenamer.renameIdentifiers(a.arg2, renamingMap)(expressionElementRenamer)
      )
  }

  private def forThreeParamFunction[A <: ThreeParamFunctionCallElement](
    constructor: (ExpressionElement, ExpressionElement, ExpressionElement) => A
  ): IdentifierLookupRenamer[A] = new IdentifierLookupRenamer[A] {
    override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
    ): A =
      constructor.apply(
        expressionElementRenamer.renameIdentifiers(a.arg1, renamingMap)(expressionElementRenamer),
        expressionElementRenamer.renameIdentifiers(a.arg2, renamingMap)(expressionElementRenamer),
        expressionElementRenamer.renameIdentifiers(a.arg3, renamingMap)(expressionElementRenamer)
      )
  }
}
