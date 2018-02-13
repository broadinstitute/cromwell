package wdl.draft3.transforms.parsing

import better.files.File
import scala.collection.JavaConverters._

import common.validation.ErrorOr.ErrorOr
import common.validation.Validation.TryValidation
import wdl.draft3.parser.WdlParser
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.FromAtoB
import wom.core.WorkflowSource

import scala.util.Try

object StringParser extends FromAtoB[FileParserInput, Ast] {
  override def convert(a: FileParserInput): ErrorOr[Ast] = Try {
    val parser = new WdlParser()
    val tokens = parser.lex(a.workflowSource, a.resource)
    val terminalMap = (tokens.asScala.toVector map {(_, a.workflowSource)}).toMap
    val syntaxErrorFormatter = WdlDraft3SyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }.toErrorOr
}

final case class FileParserInput(workflowSource: WorkflowSource, resource: String)

object FileParser extends FromAtoB[File, Ast] {
  override def convert(a: File): ErrorOr[Ast] = {
    val parserInput = FileParserInput(a.contentAsString, a.name)
    StringParser.convert(parserInput)
  }
}
