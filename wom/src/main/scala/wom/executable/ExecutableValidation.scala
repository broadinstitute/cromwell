package wom.executable

import common.validation.IOChecked
import common.validation.IOChecked._
import wom.callable.ExecutableCallable
import wom.executable.Executable.{DelayedCoercionFunction, InputParsingFunction, ResolvedExecutableInputs}
import wom.expression.IoFunctionSet
import wom.graph.Graph

private [executable] object ExecutableValidation {

  private [executable] def validateExecutable(entryPoint: ExecutableCallable,
                                              inputParsingFunction: InputParsingFunction,
                                              parseGraphInputs: (Graph, Map[String, DelayedCoercionFunction], IoFunctionSet) => IOChecked[ResolvedExecutableInputs],
                                              inputFile: Option[String],
                                              ioFunctions: IoFunctionSet): IOChecked[Executable] = for {
    parsedInputs <- inputFile.map(inputParsingFunction).map(_.toIOChecked).getOrElse(IOChecked.pure(Map.empty[String, DelayedCoercionFunction]))
    validatedInputs <- parseGraphInputs(entryPoint.graph, parsedInputs, ioFunctions)
  } yield Executable(entryPoint, validatedInputs)
}
