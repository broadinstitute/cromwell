package womtool.cmdline

import cromwell.core.path.Path

final case class PartialWomtoolCommandLineArguments(command: Option[WomtoolCommand] = None,
                                                    workflowSource: Option[Path] = None,
                                                    workflowInputs: Option[Path] = None,
                                                    displayOptionalInputs: Option[Boolean] = None,
                                                    highlightMode: Option[HighlightMode] = None,
                                                    listDependencies: Option[Boolean] = None
                                                   )

sealed trait ValidatedWomtoolCommandLine
final case class ParseCommandLine(workflowSource: Path) extends ValidatedWomtoolCommandLine
final case class ValidateCommandLine(workflowSource: Path,
                                     inputs: Option[Path],
                                     listDependencies: Boolean) extends ValidatedWomtoolCommandLine
final case class HighlightCommandLine(workflowSource: Path,
                                      highlightMode: HighlightMode) extends ValidatedWomtoolCommandLine
final case class InputsCommandLine(workflowSource: Path, showOptionals: Boolean) extends ValidatedWomtoolCommandLine
final case class WomtoolGraphCommandLine(workflowSource: Path) extends ValidatedWomtoolCommandLine
final case class WomtoolWdlUpgradeCommandLine(workflowSource: Path) extends ValidatedWomtoolCommandLine
final case class WomtoolWomGraphCommandLine(workflowSource: Path) extends ValidatedWomtoolCommandLine

sealed trait WomtoolCommand

object WomtoolCommand {
  case object Parse extends WomtoolCommand
  case object Validate extends WomtoolCommand
  case object Highlight extends WomtoolCommand
  case object Inputs extends WomtoolCommand
  case object Graph extends WomtoolCommand
  case object Upgrade extends WomtoolCommand
  case object WomGraph extends WomtoolCommand
}

sealed trait HighlightMode

object HighlightMode {
  case object HtmlHighlighting extends HighlightMode
  case object ConsoleHighlighting extends HighlightMode
  final case class UnrecognizedHighlightingMode(string: String) extends HighlightMode
}

