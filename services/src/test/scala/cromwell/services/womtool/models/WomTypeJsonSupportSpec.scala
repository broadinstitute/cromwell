package cromwell.services.womtool.models

import org.scalatest.{FlatSpec, Matchers}
import wom.types.{WomArrayType, WomStringType}
import WomTypeJsonSupport.womTypeEncoder
import io.circe.Json

class WomTypeJsonSupportSpec extends FlatSpec with Matchers {

  it should "test anything" in {
    womTypeEncoder.apply(WomArrayType(WomStringType)) shouldBe {
      Json.obj(
        (
          "typeName",
          Json.fromString("Array")
        ),
        (
          "arrayType",
          Json.obj(("typeName", Json.fromString("String")))
        )
      )
    }
  }

}
