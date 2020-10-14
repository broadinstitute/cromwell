package common.util

import common.assertion.CromwellTimeoutSpec
import common.util.TerminalUtil._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class TerminalUtilSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "TerminalUtil"

  it should "highlight" in {
    val Magenta = 35
    val result = highlight(Magenta, "hello world")
    val expected = "\u001b[38;5;35mhello world\u001b[0m"
    result should be(expected)
  }

  it should "create a markdown table" in {
    val header = Seq("col a", "col b")
    val rows = Seq(
      Seq("cell 1.a", "cell 1.b"),
      Seq("cell 2.a", "cell 2.b")
    )
    val result = mdTable(rows, header)
    val expected =
      """>|col a   |col b   |
         >|--------|--------|
         >|cell 1.a|cell 1.b|
         >|cell 2.a|cell 2.b|
         >""".stripMargin('>')
    result should be(expected)
  }

}
