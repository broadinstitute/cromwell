package wdl.draft3.transforms

import better.files.File
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ShortCircuitingFlatMap
import wdl.draft3.transforms.ast2wdlom.checkedFileElementFromFile
import wdl.draft3.transforms.wdlom2wom.fileElementToWomExecutable
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

package object parsing {
  implicit val checkedFileToAst = CheckedAtoB.fromCheck(FileParser.convert)
  implicit val checkedStringToAst = CheckedAtoB.fromCheck(StringParser.convert)

  // TODO: 2.11: the fromChecked version of this works fine in 2.12 but gets "Checked does not have a method map" in 2.11
  // Rather than fighting that, I'm going to do this inefficient thing for now and let's switch back to fromChecked in April
  implicit val checkedFileToWomExecutable = CheckedAtoB.fromErrorOr { (a: ExecutableMakerInputs[File]) =>
    for {
      fileElement <- checkedFileElementFromFile(a.from).toValidated
      inputs = ExecutableMakerInputs(fileElement, a.inputs)
      executable <- fileElementToWomExecutable(inputs).toValidated
    } yield executable
  }
}
