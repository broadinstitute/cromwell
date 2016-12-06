package cromwell.backend.validation

import cromwell.backend.validation.RuntimeAttributesDefault._
import cromwell.core.WorkflowOptions
import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import wdl4s.types._
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

class RuntimeAttributesDefaultSpec extends FlatSpec with Matchers {

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

    val coercionMap: Map[String, Set[WdlType]] = Map(
      "str" -> Set(WdlStringType),
      "bool" -> Set(WdlBooleanType),
      "number" -> Set(WdlIntegerType),
      "array" -> Set(WdlArrayType(WdlIntegerType))
    )

    val defaults = workflowOptionsDefault(workflowOptions, coercionMap)
    defaults.isSuccess shouldBe true
    defaults.get should contain theSameElementsAs Map(
      "str" -> WdlString("myString"),
      "bool" -> WdlBoolean.True,
      "number" -> WdlInteger(8),
      "array" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(2)))
    )
  }

  it should "only return default values if they're in the coercionMap" in {
    val workflowOptions = WorkflowOptions(JsObject(
      "default_runtime_attributes" -> JsObject(map))
    )

    val coercionMap: Map[String, Set[WdlType]] = Map(
      "str" -> Set(WdlStringType),
      "number" -> Set(WdlIntegerType)
    )

    val defaults = workflowOptionsDefault(workflowOptions, coercionMap)
    defaults.isSuccess shouldBe true
    defaults.get should contain theSameElementsAs Map(
      "str" -> WdlString("myString"),
      "number" -> WdlInteger(8)
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

    val coercionMap: Map[String, Set[WdlType]] = Map(
      "str" -> Set(WdlBooleanType),
      "bool" -> Set(WdlBooleanType),
      "number" -> Set(WdlIntegerType),
      "array" -> Set(WdlArrayType(WdlIntegerType))
    )

    val defaults = workflowOptionsDefault(workflowOptions, coercionMap)
    defaults.isFailure shouldBe true
    defaults.failed.get.getMessage shouldBe s"Failed to coerce default runtime options:\nCould not parse JsonValue ${map("str")} to valid WdlValue for runtime attribute str"
  }

  it should "fold default values" in {
    val rawAttributes = Map(
      "a" -> WdlString("a")
    )
    val default1 = Map(
      "a" -> WdlString("NOT a"),
      "b" -> WdlString("b")
    )
    val default2 = Map(
      "b" -> WdlString("NOT b"),
      "c" -> WdlString("c")
    )

    withDefaults(rawAttributes, List(default1, default2)) should contain theSameElementsAs Map(
      "a" -> WdlString("a"),
      "b" -> WdlString("b"),
      "c" -> WdlString("c")
    )
  }

  "noValueFoundFor" should "provide an invalidNel for missing values" in {
    import cats.syntax.validated._
    noValueFoundFor("myKey") shouldBe "Can't find an attribute value for key myKey".invalidNel
  }
}
