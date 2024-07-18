package wdl.draft3.transforms.linking.expression.files

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.types.{WomCompositeType, WomSingleFileType}
import wom.values.WomSingleFile

class FileEvaluatorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "FileEvaluator[ExpressionElement]"

  it should "find a file at the top level in WomObjects" in {
    val expressionElement: ExpressionElement = ObjectLiteral(Map("a_file" -> StringLiteral("moo.txt")))
    val structType = WomCompositeType(Map("a_file" -> WomSingleFileType))

    val evaluatedFiles = expressionElement.predictFilesNeededToEvaluate(inputs = Map.empty,
                                                                        ioFunctionSet = NoIoFunctionSet,
                                                                        coerceTo = structType
    )
    evaluatedFiles shouldBeValid Set(WomSingleFile("moo.txt"))
  }

  it should "find files in a string expression" in {
    val filePiece: StringPiece = StringLiteral("moo.txt")
    val anotherFilePiece: StringPiece = StringLiteral("boo.txt")
    val notAFilePiece: StringPiece = NewlineEscape
    val anotherNotAFilePiece: StringPiece = BackslashEscape

    val expressionElement: ExpressionElement =
      StringExpression(Seq(filePiece, anotherFilePiece, notAFilePiece, anotherNotAFilePiece))

    val evaluatedFiles = expressionElement.predictFilesNeededToEvaluate(inputs = Map.empty,
                                                                        ioFunctionSet = NoIoFunctionSet,
                                                                        coerceTo = WomSingleFileType
    )
    evaluatedFiles shouldBeValid Set(WomSingleFile("moo.txt"), WomSingleFile("boo.txt"))
  }
}
