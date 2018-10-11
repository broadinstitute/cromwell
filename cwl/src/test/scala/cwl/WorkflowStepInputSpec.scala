package cwl

import cats.data.NonEmptyList
import cwl.command.ParentName
import org.scalacheck.Prop._
import org.scalacheck.Properties
import shapeless.Coproduct
import wom.types._

object WorkflowStepInputSpec extends Properties("WorkflowStepInput") {
  implicit val pn = ParentName("_")

  val it = Coproduct[MyriadInputInnerType](CwlType.String)
  val stringType = Coproduct[MyriadInputType](it)

  val arraySchema = Coproduct[MyriadInputInnerType](InputArraySchema(items = stringType))

  val arrayStringType = WomArrayType(WomStringType)

  property("use single source type if no default value and no merge method is specified") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomStringType), None, None, hasDefault = false) ==
      Right(WomStringType)
  }

  property("use optional single source type if no default value and no merge method is specified") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomOptionalType(WomStringType)), None, None, hasDefault = false) ==
      Right(WomOptionalType(WomStringType))
  }

  property("use single source type if there's a default value and no merge method is specified") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomStringType), None, None, hasDefault = true) ==
      Right(WomStringType)
  }

  property("use unpacked single source type if there's a default value and no merge method is specified") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomOptionalType(WomStringType)), None, None, hasDefault = true) ==
      Right(WomStringType)
  }

  property("wrap single source type in an array if merge method is nested") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomStringType), Option(LinkMergeMethod.MergeNested), None, hasDefault = false) ==
      Right(WomArrayType(WomStringType))
  }

  property("find the closest common type to all sources if merge method is nested") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomStringType, "s#in2" -> WomIntegerType), Option(LinkMergeMethod.MergeNested), None, hasDefault = false) ==
      Right(WomMaybeEmptyArrayType(WomStringType))
  }

  property("validate array inner type against target type if merge method is flattened") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomStringType, "s#in2" -> WomStringType),
      Option(LinkMergeMethod.MergeFlattened),
      Option(arrayStringType), hasDefault = false
    ) == Right(WomArrayType(WomStringType))
  }

  property("validate type against target type if merge method is flattened") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomStringType, "s#in2" -> WomStringType),
      Option(LinkMergeMethod.MergeFlattened),
      Option(WomStringType), hasDefault = false
    ) == Right(WomArrayType(WomStringType))
  }

  property("fail if target type does not conform to source types if merge method is flattened") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomObjectType, "s#in2" -> WomObjectType),
      Option(LinkMergeMethod.MergeFlattened),
      Option(WomStringType), hasDefault = false
    ) == Left(NonEmptyList.one("could not verify that types Map(s#in -> WomObjectType, s#in2 -> WomObjectType) and the items type of the run's InputArraySchema WomStringType were compatible"))
  }

  property("fall back to the closest common type") = secure {
    WorkflowStepInput.determineMergeType(Map("s#in" -> WomBooleanType, "s#in2" -> WomIntegerType),
      Option(LinkMergeMethod.MergeFlattened),
      None, hasDefault = false
    ) == Right(WomStringType)
  }
}

