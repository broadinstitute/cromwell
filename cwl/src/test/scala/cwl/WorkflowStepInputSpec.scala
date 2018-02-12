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

object WorkflowStepInputSpec extends Properties("WorkflowStepInput") {
  implicit val pn = ParentName("_")

  val it = Coproduct[MyriadInputInnerType](CwlType.String)
  val stringType = Coproduct[MyriadInputType](it)

  val arraySchema = Coproduct[MyriadInputInnerType](InputArraySchema(items = stringType))

  val arrayStringType = Coproduct[MyriadInputType](arraySchema)
  val source = Some(Coproduct[InputSource]("s#in"))

  property("imply an array type if parameter is named in scatter operation") = secure {
    val wsi = WorkflowStepInput("h", source = source)
    WorkflowStepInput.determineType(wsi, Map.empty, Some(stringType), true) == Right(WomMaybeEmptyArrayType(WomStringType))
  }

  property("imply an array type if sink parameter is an array") = secure {
    val wsi = WorkflowStepInput("h")

    WorkflowStepInput.determineType(wsi, Map.empty, Some(arrayStringType), false) == Right(WomArrayType(WomStringType))
  }

  property("assert merge_nested is default LinkMergeMethod") = secure {
    val wsi = WorkflowStepInput("h")
    wsi.effectiveLinkMerge == LinkMergeMethod.MergeNested
  }

  property("assert source type is compatible with single element of 'items' type of the destination array parameter when merge_flattened option enabled") = secure {
    val miit = Coproduct[MyriadInputInnerType](InputArraySchema(items = stringType))
    val mit = Coproduct[MyriadInputType](miit)
    val inputSource =  Coproduct[InputSource](Array("_#h", "_#i"))
    val wsi = WorkflowStepInput("s#h", source = Some(inputSource), linkMerge = Some(LinkMergeMethod.MergeFlattened))

    val stringToType = Map("h" -> WomArrayType(WomStringType), "i" -> WomArrayType(WomStringType))

    WorkflowStepInput.
      determineType(wsi, stringToType, Some(mit), false) == Right(WomArrayType(WomStringType))
  }
}

