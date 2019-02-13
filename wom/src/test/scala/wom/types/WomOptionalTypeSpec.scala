package wom.types

import org.scalatest.{FlatSpecLike, Matchers}
import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor2}
import wom.values._
import wom.types.WomOptionalTypeSpecDefs._


class WomOptionalTypeSpec() extends WomCoercionSpec(goodCoercionTable, badCoercionTable, behaviorOf) with FlatSpecLike with Matchers {

  import TableDrivenPropertyChecks._

  val baseTypes = Table[WomOptionalType, WomOptionalType, WomType](
    ("optional type", "flat optional type", "base member type"),
    (WomOptionalType(WomIntegerType), WomOptionalType(WomIntegerType), WomIntegerType),
    (WomOptionalType(WomOptionalType(WomIntegerType)), WomOptionalType(WomIntegerType), WomIntegerType),
    (WomOptionalType(WomOptionalType(WomOptionalType(WomArrayType(WomOptionalType(WomIntegerType))))), WomOptionalType(WomArrayType(WomOptionalType(WomIntegerType))), WomArrayType(WomOptionalType(WomIntegerType)))
  )

  forAll(baseTypes) { (optType, flatOptType, baseType) =>
    it should s"get ${baseType.stableName} as the base type for ${optType.stableName}" in {
      optType.baseMemberType should be(baseType)
    }

    it should s"get ${flatOptType.stableName} as the flat optional type for ${optType.stableName}" in {
      optType.flatOptionalType should be(flatOptType)
    }
  }
}

private object WomOptionalTypeSpecDefs {
  import TableDrivenPropertyChecks._
  import spray.json._

  val goodCoercionTable = Table[Any, WomValue](
    ("fromValue", "toValue"),
    // Identity coercion for optionals
    (WomOptionalValue(WomInteger(75)), WomOptionalValue(WomInteger(75))),
    (WomOptionalValue(WomIntegerType, None), WomOptionalValue(WomIntegerType, None)),

    // Inner coercion defined:
    (WomOptionalValue(WomInteger(4)), WomOptionalValue(WomString("4"))),
    (WomOptionalValue(WomString("a.txt")), WomOptionalValue(WomSingleFile("a.txt"))),
    (WomOptionalValue(WomOptionalValue(WomString("a.txt"))), WomOptionalValue(WomOptionalValue(WomSingleFile("a.txt")))),
    (WomOptionalValue(WomIntegerType, None), WomOptionalValue(WomStringType, None)),

    // auto-boxing optionals
    (WomInteger(3), WomOptionalValue(WomInteger(3))),
    (WomInteger(3), WomOptionalValue(WomOptionalValue(WomInteger(3)))),
    (WomInteger(3), WomOptionalValue(WomOptionalValue(WomString("3")))), // Box and coerce
    (WomOptionalValue(WomInteger(3)), WomOptionalValue(WomOptionalValue(WomInteger(3)))),
    (WomOptionalValue(WomIntegerType, None), WomOptionalValue(WomOptionalType(WomIntegerType), None)),

    // flattening nested optionals
    (WomOptionalValue(WomOptionalValue(WomInteger(3))), WomOptionalValue(WomInteger(3))),
    (WomOptionalValue(WomOptionalValue(WomIntegerType, None)), WomOptionalValue(WomIntegerType, None)),

    // flattening and boxing and coercion all at once:
    (WomOptionalValue(WomOptionalValue(WomOptionalValue(WomInteger(3)))), WomOptionalValue(WomOptionalValue(WomString("3")))),
    (WomOptionalValue(WomOptionalValue(WomInteger(3))), WomOptionalValue(WomOptionalValue(WomOptionalValue(WomString("3"))))),
    (WomOptionalValue(WomOptionalValue(WomOptionalValue(WomIntegerType, None))), WomOptionalValue(WomOptionalType(WomStringType), None)),
    (WomOptionalValue(WomOptionalValue(WomIntegerType, None)), WomOptionalValue(WomOptionalType(WomOptionalType(WomStringType)), None)),

    // Javascript coercions:
    (JsNull, WomOptionalValue.none(WomOptionalType(WomIntegerType))),
    ("[1, 2, null, 4]".parseJson, WomArray(WomArrayType(WomOptionalType(WomIntegerType)),
      List(WomOptionalValue(WomInteger(1)), WomOptionalValue(WomInteger(2)), WomOptionalValue.none(WomIntegerType), WomOptionalValue(WomInteger(4))))),
    ("""
       |{
       |  "one": 1,
       |  "two": null,
       |  "three": 3
       |}
     """.stripMargin.parseJson,
      WomMap(WomMapType(WomStringType, WomOptionalType(WomIntegerType)), Map(
        WomString("one") -> WomOptionalValue(WomInteger(1)),
        WomString("two") -> WomOptionalValue.none(WomIntegerType),
        WomString("three") ->WomOptionalValue(WomInteger(3))
      ))
    )
  )

  val badCoercionTable: TableFor2[Any, WomType] = Table[Any, WomType](
    ("fromValue", "toType"),
    (WomOptionalValue(WomOptionalValue(WomFloat(75.0))), WomOptionalType(WomIntegerType)),
    (WomOptionalValue(WomInteger(3)), WomOptionalType(WomBooleanType)),
    (WomOptionalValue(WomOptionalValue(WomInteger(3))), WomOptionalType(WomOptionalType(WomBooleanType)))
  )

  val behaviorOf: String = "WomOptionalType"
}
