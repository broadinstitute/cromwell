package womtool.validate

import cromwell.core.path.Path
import wom.ResolvedImportRecord
import womtool.WomtoolMain.{SuccessfulTermination, Termination, UnsuccessfulTermination}
import womtool.input.WomGraphMaker

object Validate {

  def validate(main: Path, inputs: Option[Path], listDependencies: Boolean): Termination = {

    def workflowDependenciesMsg(workflowResolvedImports: Set[ResolvedImportRecord]) = {
      val msgPrefix = "\nList of Workflow dependencies is:\n"
      val dependenciesList = if (workflowResolvedImports.nonEmpty) workflowResolvedImports.map(_.importPath).mkString("\n") else "None"

      msgPrefix + dependenciesList
    }

    def validationSuccessMsg(workflowResolvedImports: Set[ResolvedImportRecord]): String = {
      val successMsg = "Success!"
      val dependenciesMsg = if (listDependencies) workflowDependenciesMsg(workflowResolvedImports) else ""
      successMsg + dependenciesMsg
    }

    if (inputs.isDefined) {
      WomGraphMaker.fromFiles(main, inputs) match {
        case Right(v) => SuccessfulTermination(validationSuccessMsg(v.resolvedImportRecords))
        case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
      }
    } else {
      WomGraphMaker.getBundle(main) match {
        case Right(b) => SuccessfulTermination(validationSuccessMsg(b.resolvedImportRecords))
        case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator))
      }
    }
  }
}
