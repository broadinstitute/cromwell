package cromwell.services.womtool.models

import org.scalatest.{FlatSpec, Matchers}
import wom.types.{WomArrayType, WomIntegerType, WomMapType, WomStringType}
import WomTypeJsonSupport.womTypeEncoder
import io.circe.Printer

class WomTypeJsonSupportSpec extends FlatSpec with Matchers {

  it should "encode Array[String] to JSON" in {
    womTypeEncoder.apply(WomArrayType(WomStringType)).pretty(Printer.spaces2.copy(colonLeft = "")) shouldBe {
      s"""{
         |  "typeName": "Array",
         |  "arrayType": {
         |    "typeName": "String"
         |  }
         |}""".stripMargin
    }
  }

  it should "encode Map[Int, Map[Map[String, Int], String]] to JSON" in {
    womTypeEncoder.apply(WomMapType(WomIntegerType, WomMapType(WomMapType(WomStringType, WomIntegerType), WomStringType))).pretty(Printer.spaces2.copy(colonLeft = "")) shouldBe {
      s"""{
         |  "typeName": "Map",
         |  "mapType": {
         |    "keyType": {
         |      "typeName": "Int"
         |    },
         |    "valueType": {
         |      "typeName": "Map",
         |      "mapType": {
         |        "keyType": {
         |          "typeName": "Map",
         |          "mapType": {
         |            "keyType": {
         |              "typeName": "String"
         |            },
         |            "valueType": {
         |              "typeName": "Int"
         |            }
         |          }
         |        },
         |        "valueType": {
         |          "typeName": "String"
         |        }
         |      }
         |    }
         |  }
         |}""".stripMargin
    }
  }

}
