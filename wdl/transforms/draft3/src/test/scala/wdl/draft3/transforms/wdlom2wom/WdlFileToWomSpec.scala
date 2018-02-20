package wdl.draft3.transforms.wdlom2wom

import cats.instances.either._
import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.ast2wdlom._

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
      (fileToAst andThen astToFileElement.map(ExecutableMakerInputs(_, List.empty, None)) andThen fileElementToWomExecutable).run(testCase) match {
        case Right(_) => // Great!
        case Left(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WOM: $formattedErrors")
      }
    }
  }
}
