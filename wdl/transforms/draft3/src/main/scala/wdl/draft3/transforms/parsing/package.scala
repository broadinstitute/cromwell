package wdl.draft3.transforms

import better.files.File
import common.transforms.CheckedAtoB
import wdl.draft3.parser.WdlParser.Ast

package object parsing {
  val fileToAst: CheckedAtoB[File, Ast] = CheckedAtoB.fromCheck(FileParser.convert)
  val stringToAst: CheckedAtoB[FileStringParserInput, Ast] = CheckedAtoB.fromCheck(StringParser.convert)
}
