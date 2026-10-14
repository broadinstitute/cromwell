package wdl.transforms.v1_2

import better.files.File
import common.transforms.CheckedAtoB
import wdl.v1_2.parser.WdlParser.Ast

package object parsing {
  val fileToAst: CheckedAtoB[File, Ast] = CheckedAtoB.fromCheck(FileParser.convert)
  val stringToAst: CheckedAtoB[FileStringParserInput, Ast] = CheckedAtoB.fromCheck(StringParser.convert)
}
