package wdl.draft3.transforms.parsing

import better.files.File
import common.Checked

import scala.collection.JavaConverters._
import common.validation.Validation.TryValidation
import wdl.draft3.parser.WdlParser
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.CheckedAtoB
import wom.core.WorkflowSource

import scala.util.Try

object StringParser {

  type FileParser = CheckedAtoB[FileParserInput, Ast]
  def instance: FileParser = CheckedAtoB(convert _)

  def convert(a: FileParserInput): Checked[Ast] = Try {
    val parser = new WdlParser()
    val tokens = parser.lex(a.workflowSource, a.resource)
    val terminalMap = (tokens.asScala.toVector map {(_, a.workflowSource)}).toMap
    val syntaxErrorFormatter = WdlDraft3SyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }.toChecked
}

final case class FileParserInput(workflowSource: WorkflowSource, resource: String)

object FileParser {

  type FileParser = CheckedAtoB[File, Ast]
  def instance: FileParser = CheckedAtoB(convert _)

  def convert(a: File): Checked[Ast] = {
    val parserInput = FileParserInput(a.contentAsString, a.name)
    StringParser.convert(parserInput)
  }
}
