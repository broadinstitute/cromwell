package wdl.draft3.transforms.wdlom2wom

import cats.instances.either._
import better.files.File
import common.transforms.CheckedAtoB
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.ast2wdlom._
import wom.executable.WomBundle

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
      val converter: CheckedAtoB[File, WomBundle] = fileToAst andThen astToFileElement.map(fe => FileElementAndImportResolvers(fe, List.empty)) andThen fileElementToWomBundle

      converter.run(testCase) match {
        case Right(_) => // Great!
        case Left(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WOM bundle: $formattedErrors")
      }
    }
  }
}
