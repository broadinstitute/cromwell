package wdl.transforms.cascades.parsing

import better.files.File
import common.Checked
import common.validation.Validation.TryValidation
import wdl.cascades.parser.WdlParser
import wdl.cascades.parser.WdlParser.Ast
import wom.core.WorkflowSource

import scala.jdk.CollectionConverters._
import scala.util.Try

object StringParser {

  def convert(a: FileStringParserInput): Checked[Ast] = Try {
    val parser = new WdlParser()
    val tokens = parser.lex(a.workflowSource, a.resource)
    val terminalMap = (tokens.asScala.toVector map { (_, a.workflowSource) }).toMap
    val syntaxErrorFormatter = WdlCascadesSyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }.toChecked
}

final case class FileStringParserInput(workflowSource: WorkflowSource, resource: String)

object FileParser {

  def convert(a: File): Checked[Ast] = {
    val parserInput = FileStringParserInput(a.contentAsString, a.name)
    StringParser.convert(parserInput)
  }
}
