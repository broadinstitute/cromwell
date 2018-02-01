package cwl

import cwl.WorkflowStep.Run
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import org.scalacheck.Properties
import org.scalacheck.Prop._
import shapeless.Coproduct
import wom.types.{WomArrayType, WomMaybeEmptyArrayType, WomStringType}

object WorkflowStepInputSpec extends Properties("WorkflowStepInput") {
  implicit val pn = ParentName("boo")

  val it = Coproduct[MyriadInputInnerType](CwlType.String)
  val expectedType = Coproduct[MyriadInputType](it)
  val source = Some(Coproduct[InputSource]("s#in"))

  property("imply an array type if named in scatter operation") = secure {
    val wsi = WorkflowStepInput("h", source = source)
    WorkflowStepInput.determineType(wsi, Map.empty, Some(expectedType), true) == Right(WomMaybeEmptyArrayType(WomStringType))
  }

  property("imply an array type if sink parameter is an array") = secure {
    val arrayType = Coproduct[MyriadInputType](Array(it))
    val wsi = WorkflowStepInput("h")

    WorkflowStepInput.determineType(wsi, Map.empty, Some(arrayType), false) == Right(WomArrayType(WomStringType))
  }

  property("assert compatible types")

}
