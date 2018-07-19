package wdl.draft3.transforms.linking.expression.files

import org.scalatest.{FlatSpec, Matchers}
import common.assertion.ErrorOrAssertions._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{ObjectLiteral, StringLiteral}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.types.{WomCompositeType, WomSingleFileType}
import wom.values.WomSingleFile
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator

class FileEvaluatorSpec extends FlatSpec with Matchers {

  behavior of "FileEvaluator[ExpressionElement]"

  it should "find a file at the top level in WomObjects" in {
    val expressionElement: ExpressionElement = ObjectLiteral(Map("a_file" -> StringLiteral("moo.txt")))
    val structType = WomCompositeType(Map("a_file" -> WomSingleFileType))

    val evaluatedFiles = expressionElement.predictFilesNeededToEvaluate(inputs = Map.empty, ioFunctionSet = NoIoFunctionSet, coerceTo = structType)
    evaluatedFiles shouldBeValid Set(WomSingleFile("moo.txt"))
  }
}
