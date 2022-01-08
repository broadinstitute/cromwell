package cromwell.backend.validation

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.core.WorkflowOptions
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._
import wom.types._
import wom.values._

class RuntimeAttributesDefaultSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "RuntimeAttributesDefaultSpec"

  val map = Map(
    "str" -> JsString("myString"),
    "bool" -> JsBoolean(true),
    "number" -> JsNumber(8),
    "array" -> JsArray(Vector(JsString("1"), JsString("2")))
  )

  it should "coerce workflow options from Json to WdlValues" in {
    val workflowOptions = WorkflowOptions(JsObject(
      "default_runtime_attributes" -> JsObject(map))
    )

    val coercionMap: Map[String, Set[WomType]] = Map(
      "str" -> Set(WomStringType),
      "bool" -> Set(WomBooleanType),
      "number" -> Set(WomIntegerType),
      "array" -> Set(WomArrayType(WomIntegerType))
    )

    val defaults = workflowOptionsDefault(workflowOptions, coercionMap)
    defaults.isSuccess shouldBe true
    defaults.get.toList should contain theSameElementsAs Map(
      "str" -> WomString("myString"),
      "bool" -> WomBoolean.True,
      "number" -> WomInteger(8),
      "array" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(2)))
    )
  }

  it should "only return default values if they're in the coercionMap" in {
    val workflowOptions = WorkflowOptions(JsObject(
      "default_runtime_attributes" -> JsObject(map))
    )

    val coercionMap: Map[String, Set[WomType]] = Map(
      "str" -> Set(WomStringType),
      "number" -> Set(WomIntegerType)
    )

    val defaults = workflowOptionsDefault(workflowOptions, coercionMap)
    defaults.isSuccess shouldBe true
    defaults.get.toList should contain theSameElementsAs Map(
      "str" -> WomString("myString"),
      "number" -> WomInteger(8)
    )
  }

  it should "return an empty map if there is no runtime attributes section in the workflow options" in {
    val workflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

    val defaults = workflowOptionsDefault(workflowOptions, Map.empty)
    defaults.isSuccess shouldBe true
    defaults.get shouldBe empty
  }

  it should "throw an exception if a value can't be coerced" in {
    val workflowOptions = WorkflowOptions(JsObject(
      "default_runtime_attributes" -> JsObject(map))
    )

    val coercionMap: Map[String, Set[WomType]] = Map(
      "str" -> Set(WomBooleanType),
      "bool" -> Set(WomBooleanType),
      "number" -> Set(WomIntegerType),
      "array" -> Set(WomArrayType(WomIntegerType))
    )

    val defaults = workflowOptionsDefault(workflowOptions, coercionMap)
    defaults.isFailure shouldBe true
    defaults.failed.get.getMessage shouldBe s"Failed to coerce default runtime options:\nCould not parse JsonValue ${map("str")} to valid WomValue for runtime attribute str"
  }

  it should "fold default values" in {
    val rawAttributes = Map(
      "a" -> WomString("a")
    )
    val default1 = Map(
      "a" -> WomString("NOT a"),
      "b" -> WomString("b")
    )
    val default2 = Map(
      "b" -> WomString("NOT b"),
      "c" -> WomString("c")
    )

    withDefaults(rawAttributes, List(default1, default2)).toList should contain theSameElementsAs Map(
      "a" -> WomString("a"),
      "b" -> WomString("b"),
      "c" -> WomString("c")
    )
  }

  "noValueFoundFor" should "provide an invalidNel for missing values" in {
    import cats.syntax.validated._
    noValueFoundFor("myKey") shouldBe "Can't find an attribute value for key myKey".invalidNel
  }
}
