package wdl.util

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class StringUtilSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "StringUtilSpec"

  val normalizeTests = Table(
    ("description", "input", "expected"),
    (
      "normalize a string with an empty line",
      "    echo hello \\\n\n    echo goodbye\n",
      "echo hello \\\n\necho goodbye"
    ),
    (
      "normalize a string with a line of spaces",
      "    echo hello \\\n    \n    echo goodbye\n",
      "echo hello \\\n\necho goodbye"
    ),
    (
      "normalize a string with extra spaces",
      "    echo hello \\\n      \n    echo goodbye\n",
      "echo hello \\\n  \necho goodbye"
    ),
    (
      "normalize a windows string with an empty line",
      "    echo hello \\\r\n\r\n    echo goodbye\n",
      "echo hello \\\n\necho goodbye"
    ),
    (
      "normalize a windows string with a line of spaces",
      "    echo hello \\\r\n    \r\n    echo goodbye\n",
      "echo hello \\\n\necho goodbye"
    ),
    (
      "normalize a windows string with extra spaces",
      "    echo hello \\\r\n      \r\n    echo goodbye\n",
      "echo hello \\\n  \necho goodbye"
    ),
    (
      "normalize a tabbed string with an empty line",
      "\t\techo hello \\\n\n\t\techo goodbye\n",
      "echo hello \\\n\necho goodbye"
    ),
    (
      "normalize a tabbed string with a line of tabs",
      "\t\techo hello \\\n\t\t\n\t\techo goodbye\n",
      "echo hello \\\n\necho goodbye"
    ),
    (
      "normalize a tabbed string with extra tabs",
      "\t\techo hello \\\n\t\t\t\n\t\techo goodbye\n",
      "echo hello \\\n\t\necho goodbye"
    ),
    (
      "normalize a string withs staggered leading spaces",
      "    echo hello \\\n  \necho goodbye\n",
      "    echo hello \\\n  \necho goodbye"
    )
  )

  forAll(normalizeTests) { (description, input, expected) =>
    it should description in {
      val actual = StringUtil.normalize(input)
      actual should be(expected)
    }
  }

  it should "not modify strings that contain only ascii characters" in {
    val input = "hi there!?"
    StringUtil.cleanUtf8mb4(input) shouldBe input
  }

  it should "not modify strings with 3-byte unicode characters" in {
    val input = "Here is my non-ascii character: \u1234 Do you like it?"
    StringUtil.cleanUtf8mb4(input) shouldBe input
  }

  it should "replace 4-byte unicode characters" in {
    val cry = new String(Character.toChars(Integer.parseInt("1F62D", 16)))
    val barf = new String(Character.toChars(Integer.parseInt("1F92E", 16)))
    val input = s"When I try to put an emoji in the database it $barf and then I $cry"
    val cleaned = "When I try to put an emoji in the database it \uFFFD and then I \uFFFD"
    StringUtil.cleanUtf8mb4(input) shouldBe cleaned
  }
}
