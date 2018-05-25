package wom.executable

import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import wom.callable.ExecutableCallable
import wom.executable.Executable.{DelayedCoercionFunction, InputParsingFunction, ResolvedExecutableInputs}
import wom.expression.IoFunctionSet
import wom.graph.Graph

private [executable] object ExecutableValidation {

  private [executable] def validateExecutable(entryPoint: ExecutableCallable,
                                              inputParsingFunction: InputParsingFunction,
                                              parseGraphInputs: (Graph, Map[String, DelayedCoercionFunction], IoFunctionSet) => ErrorOr[ResolvedExecutableInputs],
                                              inputFile: Option[String],
                                              ioFunctions: IoFunctionSet): Checked[Executable] = for {
    parsedInputs <- inputFile.map(inputParsingFunction).getOrElse(Map.empty[String, DelayedCoercionFunction].validNelCheck)
    validatedInputs <- parseGraphInputs(entryPoint.graph, parsedInputs, ioFunctions).toEither
  } yield Executable(entryPoint, validatedInputs)
}
