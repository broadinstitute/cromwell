package womtool.validate

import cromwell.core.path.Path
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

object Validate {
  def validate(main: Path, inputs: Option[Path]): Termination = {
    WomGraphMaker.fromFiles(main, inputs) match {
      case Right(_) => SuccessfulTermination("")
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
    }
  }
}
