package wom.types

import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor2}
import org.scalatest.{FlatSpecLike, Matchers}
import wom.types.WomCompositeTypeSpecDefs._
import wom.values._


class WomCompositeTypeSpec() extends WomCoercionSpec(goodCoercionTable, badCoercionTable, behaviorOf) with FlatSpecLike with Matchers {

  "WomObject with composite type" should "fail to build invalid values" in {
    val wrongType = the[Exception] thrownBy {
      WomObject.withTypeUnsafe(
        Map(
          "a" -> WomString("a")
        ),
        WomCompositeType(
          Map(
            "a" -> WomBooleanType
          )
        )
      )
    }

    wrongType.getMessage shouldBe "Error(s):\nNo coercion defined from wom value(s) '\"a\"' of type 'String' to 'Boolean'."

    val missingValues = the[Exception] thrownBy {
      WomObject.withTypeUnsafe(
        Map.empty,
        WomCompositeType(
          Map(
            "a" -> WomBooleanType
          )
        )
      )
    }

    missingValues.getMessage shouldBe "Error(s):\nNo value for field 'a' with non optional type 'Boolean' has been provided"
  }
}

private object WomCompositeTypeSpecDefs {
  import TableDrivenPropertyChecks._
  import spray.json._
  val stringIntCompositeType = WomCompositeType(Map("a" -> WomStringType, "b" -> WomIntegerType))
  val nestedStringIntCompositeType = WomCompositeType(Map("nested_a" -> WomStringType, "nested_b" -> WomIntegerType))
  val complexNestedCompositeType = WomCompositeType(Map("a" -> WomArrayType(WomStringType), "b" -> nestedStringIntCompositeType))
  val simpleCompositeValue = WomObject.withTypeUnsafe(
    Map("a" -> WomString("0"), "b" -> WomInteger(1)),
    stringIntCompositeType
  )
  val complexCompositeValue = WomObject.withTypeUnsafe(
    Map(
      // Field "a" is an array of string
      "a" -> WomArray(WomArrayType(WomStringType), List(WomString("5"))),
      // Field "b" is a composite type itself
      "b" -> WomObject.withTypeUnsafe(Map("nested_a" -> WomString("8"), "nested_b" -> WomInteger(2)), nestedStringIntCompositeType)
    ),
    complexNestedCompositeType
  )
  
  val simpleCompositeValueAsAMap = WomMap(
    WomMapType(WomStringType, WomIntegerType),
    Map(
      WomString("a") -> WomInteger(0),
      WomString("b") -> WomInteger(1)
    )
  )
  
  val goodCoercionTable = Table[Any, WomValue](
    ("fromValue", "toValue"),
    // Identity coercion with simple field types
    (simpleCompositeValue, simpleCompositeValue),
    // Identity coercion with nested complex types
    (complexCompositeValue, complexCompositeValue),
    // Coercion from WomMap
    (simpleCompositeValueAsAMap, simpleCompositeValue),
    // Coercion from untyped WomObject
    (WomObject(
      Map(
        "a" -> WomString("0"),
        "b" -> WomInteger(1)
      )
    ), simpleCompositeValue),
    (WomObject(
      Map(
        "a" -> WomArray(WomArrayType(WomIntegerType), List(WomInteger(5))),
        "b" -> WomObject(
          Map(
            "nested_a" -> WomInteger(8),
            "nested_b" -> WomInteger(2)
          )
        )
      )
    ), complexCompositeValue),
    // Coercion from Json
    ("""{ "a": "0", "b": "1" }""".parseJson, simpleCompositeValue),
    ("""{ "a": ["5"], "b": { "nested_a": "8", "nested_b": 2 } }""".parseJson, complexCompositeValue),
    // Inner values coercion (field types do not match the composite types but can be coerced to it)
    (
      WomObject.withTypeUnsafe(
        Map(
          // Array of int will be coerced to array of strings
          "a" -> WomArray(WomArrayType(WomIntegerType), List(WomInteger(5))),
          // Instead of a composite we have a Map that can be coerced to the expected composite type
          "b" -> WomMap(
            WomMapType(WomSingleFileType, WomIntegerType),
            Map(
              WomSingleFile("nested_a") -> WomInteger(8) , WomSingleFile("nested_b") -> WomInteger(2)
            )
          )
        ),
        WomCompositeType(
          Map(
            "a" -> WomArrayType(WomIntegerType),
            "b" -> WomMapType(WomSingleFileType, WomIntegerType)
          )
        )
      ),
      complexCompositeValue
    ),
    // Automatic optional none default
    (
      WomObject.withTypeUnsafe(Map.empty, WomCompositeType(Map("a" -> WomOptionalType(WomStringType)))),
      WomObject.withTypeUnsafe(Map("a" -> WomOptionalValue.none(WomStringType)), WomCompositeType(Map("a" -> WomOptionalType(WomStringType))))
    )
  )

  val badCoercionTable: TableFor2[Any, WomType] = Table[Any, WomType](
    ("fromValue", "toType"),
    // Missing fields
    (
      WomObject.withTypeUnsafe(Map.empty, WomCompositeType(Map.empty)),
      stringIntCompositeType
    ),
    // Non coerceable field types
    (
      WomObject.withTypeUnsafe(
        Map(
          "a" -> WomString("0"),
          "b" -> WomBoolean(false)
        ), WomCompositeType(Map("a" -> WomStringType, "b" -> WomBooleanType))
      ),
      stringIntCompositeType
    ),
    (
      WomMap(WomMapType(WomArrayType(WomStringType), WomStringType), Map.empty),
      stringIntCompositeType
    )
  )

  val behaviorOf: String = "WomCompositeType"
}
