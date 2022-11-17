package cwl

import cats.data.NonEmptyList
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.types._

class WorkflowStepInputSpec extends AnyFlatSpec with Matchers {
  behavior of "WorkflowStepInput"

  it should "use single source type if no default value and no merge method is specified" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomStringType),
      linkMerge = None,
      expectedTypeAsWom = None,
      hasDefault = false
    ) shouldBe
      Right(WomStringType)
  }

  it should "use optional single source type if no default value and no merge method is specified" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomOptionalType(WomStringType)),
      linkMerge = None,
      expectedTypeAsWom = None,
      hasDefault = false,
    ) shouldBe
      Right(WomOptionalType(WomStringType))
  }

  it should "use single source type if there's a default value and no merge method is specified" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomStringType),
      linkMerge = None,
      expectedTypeAsWom = None,
      hasDefault = true,
    ) shouldBe
      Right(WomStringType)
  }

  it should "use unpacked single source type if there's a default value and no merge method is specified" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomOptionalType(WomStringType)),
      linkMerge = None,
      expectedTypeAsWom = None,
      hasDefault = true,
    ) shouldBe
      Right(WomStringType)
  }

  it should "wrap single source type in an array if merge method is nested" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomStringType),
      linkMerge = Option(LinkMergeMethod.MergeNested),
      expectedTypeAsWom = None,
      hasDefault = false,
    ) shouldBe
      Right(WomArrayType(WomStringType))
  }

  it should "find the closest common type to all sources if merge method is nested" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomStringType, "s#in2" -> WomIntegerType),
      linkMerge = Option(LinkMergeMethod.MergeNested),
      expectedTypeAsWom = None,
      hasDefault = false,
    ) shouldBe
      Right(WomMaybeEmptyArrayType(WomStringType))
  }

  it should "validate array inner type against target type if merge method is flattened" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomStringType, "s#in2" -> WomStringType),
      linkMerge = Option(LinkMergeMethod.MergeFlattened),
      expectedTypeAsWom = Option(WomArrayType(WomStringType)),
      hasDefault = false,
    ) shouldBe Right(WomArrayType(WomStringType))
  }

  it should "validate type against target type if merge method is flattened" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomStringType, "s#in2" -> WomStringType),
      linkMerge = Option(LinkMergeMethod.MergeFlattened),
      expectedTypeAsWom = Option(WomStringType),
      hasDefault = false,
    ) shouldBe Right(WomArrayType(WomStringType))
  }

  it should "fail if target type does not conform to source types if merge method is flattened" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomObjectType, "s#in2" -> WomObjectType),
      linkMerge = Option(LinkMergeMethod.MergeFlattened),
      expectedTypeAsWom = Option(WomStringType),
      hasDefault = false,
    ) shouldBe Left(NonEmptyList.one(
      "could not verify that types Map(s#in -> WomObjectType, s#in2 -> WomObjectType)" +
        " and the items type of the run's InputArraySchema WomStringType were compatible"
    ))
  }

  it should "fall back to the closest common type" in {
    WorkflowStepInput.determineMergeType(
      sources = Map("s#in" -> WomBooleanType, "s#in2" -> WomIntegerType),
      linkMerge = Option(LinkMergeMethod.MergeFlattened),
      expectedTypeAsWom = None,
      hasDefault = false,
    ) shouldBe Right(WomStringType)
  }
}
