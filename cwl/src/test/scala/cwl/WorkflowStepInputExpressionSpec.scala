package cwl

import cats.data.NonEmptyList
import cats.data.Validated.Valid
import cats.syntax.either._
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import org.scalacheck.Prop._
import org.scalacheck.Properties
import shapeless.Coproduct
import wom.graph.GraphNodePort.GraphNodeOutputPort
import wom.graph.WomIdentifier
import wom.types._
import wom.values.{WomArray, WomInteger, WomString}

object WorkflowStepInputExpressionSpec extends Properties("Workflow Step Input Expression"){

  implicit val pn = ParentName("_")

  val it = Coproduct[MyriadInputInnerType](CwlType.String)

  property("assert source arrays are concatenated when merge_flattened LinkMergeMethod is used") = secure {
    val arrayType = Coproduct[MyriadInputType](Array(it))
    val inputSource =  Coproduct[InputSource](Array("i1", "i2"))
    val i1OutputPort = GraphNodeOutputPort(WomIdentifier("i1"), WomStringType, null)
    val i2OutputPort = GraphNodeOutputPort(WomIdentifier("i2"), WomStringType, null)
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource), linkMerge = Some(LinkMergeMethod.MergeFlattened))

    val mergeTpeEither = WorkflowStepInput.determineMergeType(Map("i1" -> WomArrayType(WomStringType), "i2" -> WomArrayType(WomStringType)), Some(LinkMergeMethod.MergeFlattened), Some(arrayType))
    
    val tpe: WomType = mergeTpeEither.getOrElse(throw new RuntimeException("expected a womType but evaluation of determineType failed"))

    val expression = WorkflowStepInputMergeExpression(wsi, tpe, NonEmptyList.of("i1" -> i1OutputPort, "i2" -> i2OutputPort), Vector.empty)

    val r = expression.evaluateValue(Map("i1" -> WomString("1"), "i2" -> WomString("2")), null)
    r == Valid(WomArray(Seq(WomString("1"), WomString("2"))))
  }

  property("array of one entry for each input link when merge_nested is used") = secure {
    val inputSource =  Coproduct[InputSource](Array("i1", "i2"))
    val i1OutputPort = GraphNodeOutputPort(WomIdentifier("i1"), WomStringType, null)
    val i2OutputPort = GraphNodeOutputPort(WomIdentifier("i2"), WomStringType, null)
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource))
    val expression = WorkflowStepInputMergeExpression(wsi, null, NonEmptyList.of("i1" -> i1OutputPort, "i2" -> i2OutputPort), Vector.empty)

    expression.evaluateValue(Map("i1" -> WomInteger(1), "i2" -> WomInteger(2)), null) ==
      Valid(WomArray(Seq(
        WomInteger(1),
        WomInteger(2)
      )))
  }

  property("list of one entry for when there is only one input link and when merge_nested is used") = secure {
    val inputSource =  Coproduct[InputSource]("i1")
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource))
    val i1OutputPort = GraphNodeOutputPort(WomIdentifier("i1"), WomStringType, null)
    val expression = WorkflowStepInputMergeExpression(wsi, null, NonEmptyList.of("i1" -> i1OutputPort), Vector.empty)

    expression.evaluateValue(Map("i1" -> WomInteger(1), "i2" -> WomInteger(2)), null) ==
      Valid(WomInteger(1))
  }

}
