package cwl

import cats.syntax.validated._
import cats.syntax.either._
import io.circe.literal._
import io.circe.{Json, yaml}
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import wom.callable.{CallableTaskDefinition, ExecutableCallable, ExecutableTaskDefinition}
import wom.executable.Executable
import wom.executable.Executable.{InputParsingFunction, ParsedInputMap}

// WARNING! Because of 2.11 vs 2.12 incompatibilities, there are two versions of this file.
// If you're making changes here, you'll also need to update ../../scala-2.12/cwl/CwlExecutableValidation.scala
// (ExecutableValidation.scala has more info on why this was necessary)
object CwlExecutableValidation {

  // Decodes the input file, and build the ParsedInputMap
  private val inputCoercionFunction: InputParsingFunction =
    inputFile => {
      yaml.parser.parse(inputFile).flatMap(_.as[Map[String, Json]]) match {
        case Left(error) => error.getMessage.invalidNelCheck[ParsedInputMap]
        case Right(inputValue) => inputValue.map({ case (key, value) => key -> value.foldWith(CwlJsonToDelayedCoercionFunction) }).validNelCheck
      }
    }

  def buildWomExecutable(callable: Checked[ExecutableCallable], inputFile: Option[String]): Checked[Executable] = {
    for {
      womDefinition <- callable
      executable <- Executable.withInputs(womDefinition, inputCoercionFunction, inputFile)
    } yield executable
  }

  def buildWomExecutable(callableTaskDefinition: ErrorOr[CallableTaskDefinition], inputFile: Option[String]): Checked[Executable] = {
    import cats.data.Validated._
    import common.validation.ErrorOr._
    (for {
      taskDefinition <- callableTaskDefinition
      executableTaskDefinition = ExecutableTaskDefinition.tryApply(taskDefinition)
      executable <- fromEither(CwlExecutableValidation.buildWomExecutable(executableTaskDefinition.toEither, inputFile))
    } yield executable).toEither
  }
}
