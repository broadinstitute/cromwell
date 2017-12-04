package cromwell.backend

import org.scalacheck.{Arbitrary, Gen, Properties}
import org.scalacheck.Prop._
import wom.expression.WomExpression
import wom.values.{WomString, WomValue}

class RuntimeAttributeValidationSpec extends Properties("Runtime Attribute Validation") {

  import WomGenerators._

  property("use default and validate it when runtime is not specified") = forAll {
    (taskName: String, attributeName: String, womValue: WomValue) =>
      val defaultRuntimeAttributes = Map(attributeName -> womValue)

      val defaultValue = womValue.asWomExpression
      val validator: Option[WomExpression] => Boolean = _.contains(defaultValue)
      BackendWorkflowInitializationActor.validateRuntimeAttributes(taskName, defaultRuntimeAttributes, Map.empty, Map((attributeName,validator)) ).isValid
  }

  property("return invalid if validator fails the test ") = forAll {
    (taskName: String, attributeName: String, womValue: WomValue) =>
      val defaultRuntimeAttributes = Map(attributeName -> womValue)

      BackendWorkflowInitializationActor.validateRuntimeAttributes(taskName, defaultRuntimeAttributes, Map.empty, Map((attributeName,(_: Option[WomExpression]) => false))).isInvalid
  }

  property("use runtime setting (not default) when both are set") = forAll {
    (taskName: String, attributeName: String, defaultWomValue: WomValue, runtimeWomExpression: WomExpression) =>
      val defaultRuntimeAttributes = Map(attributeName -> defaultWomValue)
      val runtimeAttributes = Map(attributeName -> runtimeWomExpression)

      val validator: Option[WomExpression] => Boolean = _.contains(runtimeWomExpression)
      BackendWorkflowInitializationActor.validateRuntimeAttributes(taskName, defaultRuntimeAttributes, runtimeAttributes, Map((attributeName,validator))).isValid
  }

  property("fail validation if no setting is present but it should be") = forAll {
    (taskName: String, attributeName: String) =>

      val validator: Option[WomExpression] => Boolean = {
        case None => false
        case Some(x) => throw new RuntimeException(s"expecting the runtime validator to receive a None but got $x")
      }
      BackendWorkflowInitializationActor.validateRuntimeAttributes(taskName, Map.empty, Map.empty, Map((attributeName,validator))).isInvalid
  }

  property("use the taskName and attribute name in correct places for failures") = forAll {
    (taskName: String, attributeName: String) =>

      val validator: Option[WomExpression] => Boolean = {
        case None => false
        case Some(x) => throw new RuntimeException(s"expecting the runtime validator to receive a None but got $x")
      }
      BackendWorkflowInitializationActor.validateRuntimeAttributes(taskName, Map.empty, Map.empty, Map((attributeName,validator))).fold(
        { errors =>
          val error = errors.toList.head
          all(
            "attribute name should be set correctly" |: error.runtimeAttributeName == attributeName,
            "task name should be set correctly" |: error.jobTag == taskName
          )
        },
        _ =>  "expecting validation to fail!" |: false
      )
  }
}

object WomGenerators {
  implicit def womValue: Arbitrary[WomValue] = Arbitrary(Gen.oneOf(WomString("some"), WomString("other")))

  implicit def womExpression: Arbitrary[WomExpression] = Arbitrary(womValue.arbitrary.map(_.asWomExpression))
}
