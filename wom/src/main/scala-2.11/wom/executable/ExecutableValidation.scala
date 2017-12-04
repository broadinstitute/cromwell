package wom.executable

import cats.syntax.either._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import wom.callable.ExecutableCallable
import wom.executable.Executable.{DelayedCoercionFunction, InputParsingFunction, ResolvedExecutableInputs}
import wom.graph.Graph

// WARNING! Because of 2.11 vs 2.12 incompatibilities, there are two versions of this file.
// If you're making changes here, you'll also need to update ../../../scala-2.12/wom/executable/ExecutableValidation.scala
// (More info in the 2.12 version on why this was necessary)
private [executable] object ExecutableValidation {

  private [executable] def validateExecutable(entryPoint: ExecutableCallable,
                                              inputParsingFunction: InputParsingFunction,
                                              parseGraphInputs: (Graph, Map[String, DelayedCoercionFunction]) => ErrorOr[ResolvedExecutableInputs],
                                              inputFile: Option[String]): Checked[Executable] = for {
    parsedInputs <- inputFile.map(inputParsingFunction).getOrElse(Map.empty[String, DelayedCoercionFunction].validNelCheck)
    validatedInputs <- parseGraphInputs(entryPoint.graph, parsedInputs).toEither
  } yield Executable(entryPoint, validatedInputs)
}
