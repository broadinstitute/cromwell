package wdl.draft3.transforms.wdlom2wdl

import java.nio.file.Paths

import cats.implicits._
import common.Checked
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.transforms.ast2wdlom.astToFileElement
import wdl.draft3.transforms.parsing.{FileStringParserInput, fileToAst, stringToAst}
import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
import wdl.model.draft3.elements.FileElement

class WdlomToWdlFileSpec extends FlatSpec with Matchers {

  it should "write out a file that re-evaluates into the same case class structure" in {
    val file = Paths.get("wdl/transforms/draft3/src/test/cases/simple_first_test.wdl")

    val model: Checked[FileElement] = (fileToAst andThen astToFileElement).run(file)

    model match {
      case Right(wdlModel) =>

        val newModel = (stringToAst andThen astToFileElement).run(FileStringParserInput(wdlModel.toWdl, "simple_first_test.wdl"))

        // Scala case class deep equality is so nice here
        newModel shouldEqual model
      case Left(_) => fail("Could not load original AST")
    }
  }

}
