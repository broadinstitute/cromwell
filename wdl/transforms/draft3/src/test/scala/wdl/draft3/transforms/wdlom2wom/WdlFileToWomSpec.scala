package wdl.draft3.transforms.wdlom2wom

import better.files.File
import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}

class WdlFileToWomSpec extends FlatSpec with Matchers {
  behavior of "WDL File to WOM"

  val testCases = File("wdl/transforms/draft3/src/test/cases")

  it should "be set up for testing" in {
    testCases.exists shouldBe true
    testCases.list.nonEmpty shouldBe true
  }

  testCases.list.filter(x => x.isRegularFile && x.extension.contains(".wdl")) foreach { testCase =>

    val fileName = testCase.name

    val itShouldString = s"create a valid WOM object for $fileName"
    val testOrIgnore: (=>Any) => Unit = if (testCase.name.endsWith(".ignored.wdl") || testCase.name.endsWith(".nowom.wdl")) {
      (it should itShouldString).ignore _
    } else {
      (it should itShouldString).in _
    }

    testOrIgnore {

      womFromDraft3FileConverter(None).convert(testCase) match {
        case Valid(_) => // Great!
        case Invalid(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WOM: $formattedErrors")
      }
    }
  }
}
