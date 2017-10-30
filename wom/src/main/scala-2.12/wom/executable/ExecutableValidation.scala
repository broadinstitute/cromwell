package wom.executable

import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.Checked._
import wom.callable.ExecutableCallable
import wom.executable.Executable.{DelayedCoercionFunction, InputParsingFunction, ResolvedExecutableInputs}
import wom.graph.Graph

// WARNING! Because of 2.11 vs 2.12 incompatibilities, there are two versions of this file.
// If you're making changes here, you'll also need to update ../../../scala-2.11/wom/executable/ExecutableValidation.scala
// (More info below on why this was necessary)
private [executable] object ExecutableValidation {
  /*
    * Validates and build an executable.
    *
    * There is a 2.12 and 2.11 version of this function because scala 2.12 has right biased Either whereas
    * scala 2.11 has not. Because of this we can map and flatMap Eithers in scala 2.12 without any external library,
    * but we need cats.syntax.either._ in 2.11. Because we warn on unused imports and warnings are fatal,
    * the same implementation cannot compile on 2.11 and 2.12 at the same time.
   */
  private [executable] def validateExecutable(entryPoint: ExecutableCallable,
                                              inputParsingFunction: InputParsingFunction,
                                              parseGraphInputs: (Graph, Map[String, DelayedCoercionFunction]) => ErrorOr[ResolvedExecutableInputs],
                                              inputFile: Option[String]): Checked[Executable] = for {
    parsedInputs <- inputFile.map(inputParsingFunction).getOrElse(Map.empty[String, DelayedCoercionFunction].validNelCheck)
    validatedInputs <- parseGraphInputs(entryPoint.graph, parsedInputs).toEither
  } yield Executable(entryPoint, validatedInputs)
}
