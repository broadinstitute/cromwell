package common.numeric

import common.numeric.IntegerUtil.IntEnhanced
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class IntegerUtilSpec extends AnyFlatSpec with Matchers {
  it should "return ordinal String for any Int" in {
    val numbers = List(0, 1, 2, 3, 4,
        10, 11, 12, 13, 14,
        20, 21, 22, 23, 24,
        100, 101, 102, 103, 104) map { _.toOrdinal }

    val expected = List("0th", "1st", "2nd", "3rd", "4th",
    "10th", "11th", "12th", "13th", "14th",
    "20th", "21st", "22nd", "23rd", "24th",
    "100th", "101st", "102nd", "103rd", "104th")

    numbers should contain theSameElementsInOrderAs expected
  }

}
