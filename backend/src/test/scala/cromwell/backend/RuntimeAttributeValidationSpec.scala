package cromwell.backend

import org.scalacheck.{Arbitrary, Gen}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.scalacheck.ScalaCheckDrivenPropertyChecks
import wom.expression.WomExpression
import wom.values.{WomString, WomValue}

class RuntimeAttributeValidationSpec extends AnyFlatSpec with Matchers with ScalaCheckDrivenPropertyChecks {
  behavior of "RuntimeAttributeValidation"

  import WomGenerators._

  it should "use default and validate it when runtime is not specified" in forAll {
    (taskName: String, attributeName: String, womValue: WomValue) =>
      val defaultRuntimeAttributes = Map(attributeName -> womValue)

      val defaultValue = womValue.asWomExpression
      val validator: Option[WomExpression] => Boolean = _.contains(defaultValue)
      assert(
        BackendWorkflowInitializationActor
          .validateRuntimeAttributes(
            taskName = taskName,
            defaultRuntimeAttributes = defaultRuntimeAttributes,
            runtimeAttributes = Map.empty,
            runtimeAttributeValidators = Map((attributeName, validator))
          )
          .isValid
      )
  }

  it should "return invalid if validator fails the test" in forAll {
    (taskName: String, attributeName: String, womValue: WomValue) =>
      val defaultRuntimeAttributes = Map(attributeName -> womValue)

      assert(
        BackendWorkflowInitializationActor
          .validateRuntimeAttributes(
            taskName = taskName,
            defaultRuntimeAttributes = defaultRuntimeAttributes,
            runtimeAttributes = Map.empty,
            runtimeAttributeValidators = Map((attributeName, (_: Option[WomExpression]) => false))
          )
          .isInvalid
      )
  }

  it should "use runtime setting (not default) when both are set" in forAll {
    (taskName: String, attributeName: String, defaultWomValue: WomValue, runtimeWomExpression: WomExpression) =>
      val defaultRuntimeAttributes = Map(attributeName -> defaultWomValue)
      val runtimeAttributes = Map(attributeName -> runtimeWomExpression)

      val validator: Option[WomExpression] => Boolean = _.contains(runtimeWomExpression)
      assert(
        BackendWorkflowInitializationActor
          .validateRuntimeAttributes(
            taskName = taskName,
            defaultRuntimeAttributes = defaultRuntimeAttributes,
            runtimeAttributes = runtimeAttributes,
            runtimeAttributeValidators = Map((attributeName, validator))
          )
          .isValid
      )
  }

  it should "fail validation if no setting is present but it should be" in forAll {
    (taskName: String, attributeName: String) =>
      val validator: Option[WomExpression] => Boolean = {
        case None => false
        case Some(x) => throw new RuntimeException(s"expecting the runtime validator to receive a None but got $x")
      }
      assert(
        BackendWorkflowInitializationActor
          .validateRuntimeAttributes(
            taskName = taskName,
            defaultRuntimeAttributes = Map.empty,
            runtimeAttributes = Map.empty,
            runtimeAttributeValidators = Map((attributeName, validator))
          )
          .isInvalid
      )
  }

  it should "use the taskName and attribute name in correct places for failures" in forAll {
    (taskName: String, attributeName: String) =>
      val validator: Option[WomExpression] => Boolean = {
        case None => false
        case Some(x) => throw new RuntimeException(s"expecting the runtime validator to receive a None but got $x")
      }
      BackendWorkflowInitializationActor
        .validateRuntimeAttributes(taskName, Map.empty, Map.empty, Map((attributeName, validator)))
        .fold(
          { errors =>
            val error = errors.toList.head
            withClue("attribute name should be set correctly")(error.runtimeAttributeName shouldBe attributeName)
            withClue("task name should be set correctly")(error.jobTag shouldBe taskName)
          },
          _ => fail("expecting validation to fail!")
        )
  }
}

object WomGenerators {
  implicit def womValue: Arbitrary[WomValue] = Arbitrary(Gen.oneOf(WomString("some"), WomString("other")))

  implicit def womExpression: Arbitrary[WomExpression] = Arbitrary(womValue.arbitrary.map(_.asWomExpression))
}
