package wdl.transforms.cascades

import better.files.File
import common.transforms.CheckedAtoB
import wdl.cascades.parser.WdlParser.Ast

package object parsing {
  val fileToAst: CheckedAtoB[File, Ast] = CheckedAtoB.fromCheck(FileParser.convert)
  val stringToAst: CheckedAtoB[FileStringParserInput, Ast] = CheckedAtoB.fromCheck(StringParser.convert)
}
