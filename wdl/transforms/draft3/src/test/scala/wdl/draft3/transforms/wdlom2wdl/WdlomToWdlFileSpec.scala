package wdl.draft3.transforms.wdlom2wdl

import better.files.File
import cats.implicits._
import common.Checked
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.transforms.ast2wdlom.astToFileElement
import wdl.draft3.transforms.parsing.{FileStringParserInput, fileToAst, stringToAst}
import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
import wdl.model.draft3.elements.FileElement

class WdlomToWdlFileSpec extends FlatSpec with Matchers {

  val testCases = File("wdl/transforms/draft3/src/test/cases")

  testCases.list.foreach { file =>

    it should s"write a file that re-evaluates to the same case classes for ${file.name}" in {

      val model: Checked[FileElement] = (fileToAst andThen astToFileElement).run(file)
        model match {
        case Right(wdlModel) =>

          val newModel = (stringToAst andThen astToFileElement).run(FileStringParserInput(wdlModel.toWdl, file.name))

          // Scala case class deep equality is so nice here
          newModel shouldEqual model
        case Left(_) => fail("Could not load original")
      }
    }
  }
}
