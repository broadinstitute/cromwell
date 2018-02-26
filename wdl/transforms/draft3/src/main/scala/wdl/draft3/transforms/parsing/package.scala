package wdl.draft3.transforms

import better.files.File
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ShortCircuitingFlatMap
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.fileToFileElement
import wdl.draft3.transforms.wdlom2wom.fileElementToWomExecutable
import wom.executable.Executable
import wom.transforms.WomExecutableMaker.WomBundleMakerInputs

package object parsing {
  val fileToAst: CheckedAtoB[File, Ast] = CheckedAtoB.fromCheck(FileParser.convert)
  val stringToAst: CheckedAtoB[FileStringParserInput, Ast] = CheckedAtoB.fromCheck(StringParser.convert)
}
