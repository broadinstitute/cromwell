package cwl

import cats.syntax.validated._
import cats.syntax.either._
import io.circe._
import io.circe.shapes._
import io.circe.generic.auto._
import eu.timepit.refined.string._
import io.circe.refined._
import io.circe.yaml
import io.circe.literal._
import lenthall.Checked
import lenthall.validation.Checked._
import wom.callable.Callable
import wom.executable.Executable
import wom.executable.Executable.{InputParsingFunction, ParsedInputMap}

// See explanation as to why there are 2 versions of this in ExecutableValidation
object CwlExecutableValidation {

  implicit val f = implicitly[Decoder[File]]

  // Decodes the input file, and build the ParsedInputMap
  private val inputCoercionFunction: InputParsingFunction =
    inputFile => {
      yaml.parser.parse(inputFile).flatMap(_.as[Map[String, MyriadInputValue]]) match {
        case Left(error) => error.getMessage.invalidNelCheck[ParsedInputMap]
        case Right(inputValue) => inputValue.map({ case (key, value) => key -> value.fold(CwlInputCoercion) }).validNelCheck
      }
    }

  def builWomExecutable(callable: Checked[Callable], inputFile: Option[String]): Checked[Executable] = {
    for {
      womDefinition <- callable
      executable <- Executable.withInputs(womDefinition, inputCoercionFunction, inputFile)
    } yield executable
  }
}
