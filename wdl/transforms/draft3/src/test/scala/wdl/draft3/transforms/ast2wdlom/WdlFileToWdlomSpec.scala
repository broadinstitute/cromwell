package wdl.draft3.transforms.ast2wdlom

import better.files.File
import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.parser.WdlParser
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.parsing.WdlDraft3SyntaxErrorFormatter
import wdl.model.draft3.elements.{FileElement, WorkflowDefinitionElement}
import wdl.draft3.transforms.ast2wdlom.WdlFileToWdlomSpec._
import wom.core.WorkflowSource

import scala.collection.JavaConverters._

class WdlFileToWdlomSpec extends FlatSpec with Matchers {

  behavior of "WDL File to WDLOM"

  val testCases = File("wdl/transforms/draft3/src/test/cases")

  it should "be set up for testing" in {
    testCases.exists shouldBe true
    testCases.list.nonEmpty shouldBe true
  }

  testCases.list.filter(x => x.isRegularFile && x.extension.contains(".wdl")) foreach { testCase =>

    val fileName = testCase.name
    val testName = testCase.nameWithoutExtension

    val itShouldString = s"create the correct Element structure for $fileName"
    val testOrIgnore: (=>Any) => Unit = if (testCase.name.endsWith(".ignored.wdl")) {
      (it should itShouldString).ignore _
    } else {
      (it should itShouldString).in _
    }

    testOrIgnore {
      val wdlFile = testCase

      val fileContents = wdlFile.contentAsString
      val ast = parseFile(fileContents, testCase.name)

      val expected = expectations.getOrElse(testName, fail(s"No Element expectation defined for $testName"))
      FromAstNode[FileElement](ast) match {
        case Valid(actual) => actual shouldBe expected
        case Invalid(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WDLOM:$formattedErrors")
      }

    }
  }
}

object WdlFileToWdlomSpec {

  val expectations: Map[String, FileElement] = Map(
    "empty_workflow" ->
      FileElement(
        imports = List.empty,
        workflows = List(WorkflowDefinitionElement("empty")),
        tasks = List.empty),
    "passthrough_workflow" ->
      FileElement(
        imports = List.empty,
        workflows = List.empty,
        tasks = List.empty),
    "simpleFirstTest" ->
      FileElement(
        imports = List.empty,
        workflows = List.empty,
        tasks = List.empty)
  )


  private def parseFile(workflowSource: WorkflowSource, resource: String): Ast = {
    val parser = new WdlParser()
    val tokens = parser.lex(workflowSource, resource)
    val terminalMap = (tokens.asScala.toVector map {(_, workflowSource)}).toMap
    val syntaxErrorFormatter = WdlDraft3SyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }
}
