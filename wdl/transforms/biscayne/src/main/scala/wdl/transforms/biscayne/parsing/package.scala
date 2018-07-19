package wdl.transforms.biscayne

import better.files.File
import common.transforms.CheckedAtoB
import wdl.biscayne.parser.WdlParser.Ast

package object parsing {
  val fileToAst: CheckedAtoB[File, Ast] = CheckedAtoB.fromCheck(FileParser.convert)
  val stringToAst: CheckedAtoB[FileStringParserInput, Ast] = CheckedAtoB.fromCheck(StringParser.convert)
}
