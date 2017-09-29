package wdl.util

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

class StringUtilSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

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

}
