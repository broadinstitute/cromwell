package wom.executable

import cats.syntax.either._
import lenthall.Checked
import lenthall.validation.Checked._
import lenthall.validation.ErrorOr.ErrorOr
import wom.callable.Callable
import wom.executable.Executable.{DelayedCoercionFunction, InputParsingFunction, ResolvedExecutableInputs}
import wom.graph.Graph

private [executable] object ExecutableValidation {
  /*
    * Validates and build an executable.
    * Note: There is a 2.12 and 2.11 version of this function.
    * This is because scala 2.12 has right biased Either whereas scala 2.11 has not
    * Because of this we can map and flatMap Eithers in scala 2.12 without any external librabry,
    * but we need cats.syntax.either._ in 2.11. Because we warn on unused imports and warnings are fatal,
    * the same implementation cannot compile on 2.11 and 2.12 at the same time.
   */
  private [executable] def validateExecutable(entryPoint: Callable,
                                              inputParsingFunction: InputParsingFunction,
                                              parseGraphInputs: (Graph, Map[String, DelayedCoercionFunction]) => ErrorOr[ResolvedExecutableInputs],
                                              inputFile: Option[String]): Checked[Executable] = for {
    validGraph <- entryPoint.graph.toEither
    parsedInputs <- inputFile.map(inputParsingFunction).getOrElse(Map.empty[String, DelayedCoercionFunction].validNelCheck)
    validatedInputs <- parseGraphInputs(validGraph, parsedInputs).toEither
  } yield Executable(entryPoint, validatedInputs)
}
