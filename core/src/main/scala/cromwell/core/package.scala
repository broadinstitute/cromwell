package cromwell

import cats.data.Validated._
import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.values.WdlValue

import scala.util.{Failure, Success, Try}
import scala.language.implicitConversions

package object core {
  type LocallyQualifiedName = String
  type FullyQualifiedName = String
  type WorkflowOutputs = Map[FullyQualifiedName, JobOutput]
  type WorkflowOptionsJson = String
  type CallOutputs = Map[LocallyQualifiedName, JobOutput]
  type HostInputs = Map[String, WdlValue]
  type EvaluatedRuntimeAttributes = Map[String, WdlValue]

  implicit def tryToErrorOr[A](trySomething: Try[A]): ErrorOr[A] = trySomething match {
    case Success(options) => options.validNel
    case Failure(err) => err.getMessage.invalidNel
  }

  implicit def errorOrToTry[A](validatedString: ErrorOr[A]): Try[A] = validatedString match {
    case Valid(options) => Success(options)
    case Invalid(err) => Failure(new RuntimeException(s"Error(s): ${err.toList.mkString(",")}"))
  }
}
