package cwl

import cats.syntax.either._
import common.Checked
import common.validation.Checked._
import io.circe.Json
import io.circe.yaml
import wom.callable.{ExecutableCallable, TaskDefinition}
import wom.executable.Executable
import wom.executable.Executable.{InputParsingFunction, ParsedInputMap}
import wom.expression.IoFunctionSet

object CwlExecutableValidation {

  // Decodes the input file, and build the ParsedInputMap
  private val inputCoercionFunction: InputParsingFunction =
    inputFile => {
      yaml.parser.parse(inputFile).flatMap(_.as[Map[String, Json]]) match {
        case Left(error) => error.getMessage.invalidNelCheck[ParsedInputMap]
        case Right(inputValue) => inputValue.map({ case (key, value) => key -> value.foldWith(CwlJsonToDelayedCoercionFunction) }).validNelCheck
      }
    }

  def buildWomExecutableCallable(callable: Checked[ExecutableCallable], inputFile: Option[String], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = {
    for {
      womDefinition <- callable
      executable <- Executable.withInputs(womDefinition, inputCoercionFunction, inputFile, ioFunctions, strictValidation)
    } yield executable
  }

  def buildWomExecutable(callableTaskDefinition: Checked[TaskDefinition], inputFile: Option[String], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = {
    for {
      taskDefinition <- callableTaskDefinition
      executableTaskDefinition = taskDefinition.toExecutable.toEither
      executable <- CwlExecutableValidation.buildWomExecutableCallable(executableTaskDefinition, inputFile, ioFunctions, strictValidation)
    } yield executable
  }
}
