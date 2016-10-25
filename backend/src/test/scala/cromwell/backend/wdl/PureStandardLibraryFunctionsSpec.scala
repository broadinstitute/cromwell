package cromwell.backend.wdl

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.expression.PureStandardLibraryFunctions
import wdl4s.types.{WdlArrayType, WdlIntegerType}
import wdl4s.values.{WdlArray, WdlInteger}

import scala.util.Success


class PureStandardLibraryFunctionsSpec extends FlatSpec with Matchers {

  behavior of "transpose"

  it should "transpose a 2x3 into a 3x2" in {
    val inArray = WdlArray(WdlArrayType(WdlArrayType(WdlIntegerType)), List(
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(2), WdlInteger(3))),
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(4), WdlInteger(5), WdlInteger(6)))
    ))

    val expectedResult = WdlArray(WdlArrayType(WdlArrayType(WdlIntegerType)), List(
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(4))),
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(2), WdlInteger(5))),
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(3), WdlInteger(6)))
    ))

    PureStandardLibraryFunctions.transpose(Seq(Success(inArray))) should be(Success(expectedResult))
  }

}
