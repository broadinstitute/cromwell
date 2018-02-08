package cwl

import cats.data.Validated.Valid
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStep.Run
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import org.scalacheck.Properties
import org.scalacheck.Prop._
import org.scalacheck.Gen._
import shapeless.Coproduct
import wom.types._
import wom.values.{WomArray, WomInteger, WomMap, WomString, WomValue}
import cats.syntax.validated._
import cats.syntax.either._

object WorkflowStepInputExpressionSpec extends Properties("Workflow Step Input Expression"){

  implicit val pn = ParentName("_")

  val it = Coproduct[MyriadInputInnerType](CwlType.String)

  property("assert source arrays are concatenated when merge_flattened LinkMergeMethod is used") = secure {
    val arrayType = Coproduct[MyriadInputType](Array(it))
    val inputSource =  Coproduct[InputSource](Array("_#i1", "_#i2"))
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource), linkMerge = Some(LinkMergeMethod.MergeFlattened))

    val tpeEither = WorkflowStepInput.determineType(wsi, Map("i1" -> WomArrayType(WomStringType), "i2" -> WomArrayType(WomStringType)), Some(arrayType), false)

    val tpe: WomType = tpeEither.getOrElse(throw new RuntimeException("expected a womType but evaluation of determineType failed"))

    val expression = WorkflowStepInputExpression(wsi, tpe, Set.empty, Vector.empty)

    expression.evaluateValue(Map("i1" -> WomString("1"), "i2" -> WomString("2")), null) ==
      Valid(WomArray(Seq(WomString("1"), WomString("2"))))
  }

  property("array of one entry for each input link when merge_nested is used") = secure {
    val inputSource =  Coproduct[InputSource](Array("_#i1", "_#i2"))
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource))
    val expression = WorkflowStepInputExpression(wsi, null, Set.empty, Vector.empty)

    expression.evaluateValue(Map("i1" -> WomInteger(1), "i2" -> WomInteger(2)), null) ==
      Valid(WomArray(Seq(
        WomMap(WomMapType(WomStringType, WomIntegerType) ,Map(WomString("i1") -> WomInteger(1))),
        WomMap(WomMapType(WomStringType, WomIntegerType) ,Map(WomString("i2") -> WomInteger(2)))
      )))
  }

  property("list of one entry for when there is only one input link and when merge_nested is used") = secure {
    val inputSource =  Coproduct[InputSource]("_#i1")
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource))
    val expression = WorkflowStepInputExpression(wsi, null, Set.empty, Vector.empty)

    expression.evaluateValue(Map("i1" -> WomInteger(1), "i2" -> WomInteger(2)), null) ==
      Valid(WomInteger(1))
  }

}
