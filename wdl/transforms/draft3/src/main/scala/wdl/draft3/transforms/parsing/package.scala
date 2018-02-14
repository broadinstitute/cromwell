package wdl.draft3.transforms

import better.files.File
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ShortCircuitingFlatMap
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.fileToFileElement
import wdl.draft3.transforms.wdlom2wom.elementToWomExecutable
import wom.executable.Executable
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

package object parsing {
  val fileToAst: CheckedAtoB[File, Ast] = CheckedAtoB.fromCheck(FileParser.convert)
  val stringToAst: CheckedAtoB[FileStringParserInput, Ast] = CheckedAtoB.fromCheck(StringParser.convert)

  // TODO: 2.11: the fromChecked version of this works fine in 2.12 but gets "Checked does not have a method map" in 2.11
  // Rather than fighting that, I'm going to do this inefficient thing for now and let's switch back to fromChecked in April
  val fileToWomExecutable: CheckedAtoB[ExecutableMakerInputs[File], Executable] = CheckedAtoB.fromErrorOr { a =>
    for {
      fileElement <- fileToFileElement(a.from).toValidated
      inputs = ExecutableMakerInputs(fileElement, a.inputs)
      executable <- elementToWomExecutable(inputs).toValidated
    } yield executable
  }
}
