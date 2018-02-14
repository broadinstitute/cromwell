package wdl.draft3.transforms

import better.files.File
import common.transforms.CheckedAtoB
import wdl.draft3.transforms.ast2wdlom.checkedFileElementFromFile
import wdl.draft3.transforms.wdlom2wom.fileElementToWomExecutable
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

package object parsing {
  implicit val checkedFileToAst = CheckedAtoB.fromCheck(FileParser.convert)
  implicit val checkedStringToAst = CheckedAtoB.fromCheck(StringParser.convert)

  implicit val checkedFileToWomExecutable = CheckedAtoB.fromCheck { (a: ExecutableMakerInputs[File]) =>
    for {
      fileElement <- checkedFileElementFromFile(a.from)
      inputs = ExecutableMakerInputs(fileElement, a.inputs)
      executable <- fileElementToWomExecutable(inputs)
    } yield executable
  }
}
