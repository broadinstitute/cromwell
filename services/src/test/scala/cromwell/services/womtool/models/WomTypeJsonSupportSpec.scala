package cromwell.services.womtool.models

import org.scalatest.{FlatSpec, Matchers}
import wom.types.{WomArrayType, WomStringType}
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

}
