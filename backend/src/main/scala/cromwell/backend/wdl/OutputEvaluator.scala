package cromwell.backend.wdl

import cats.data.EitherT._
import cats.data.Validated.{Invalid, Valid}
import cats.data.{EitherT, NonEmptyList}
import cats.instances.try_._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.{CallOutputs, JobOutput}
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import wdl.types.WdlType
import wdl.values.WdlValue
import wom.callable.Callable.OutputDefinition
import wom.expression.IoFunctionSet

import scala.util.{Failure, Success, Try}
object OutputEvaluator {
  sealed trait EvaluatedJobOutputs
  case class ValidJobOutputs(outputs: CallOutputs) extends EvaluatedJobOutputs
  case class InvalidJobOutputs(errors: NonEmptyList[String]) extends EvaluatedJobOutputs
  case class JobOutputsEvaluationException(exception: Throwable) extends EvaluatedJobOutputs

  type OutputResult[A] = EitherT[Try, NonEmptyList[String], A]

  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      ioFunctions: IoFunctionSet,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): EvaluatedJobOutputs = {
    val knownValues = jobDescriptor.localInputs
    
    def foldFunction(accumulatedOutputs: Try[ErrorOr[List[(String, WdlValue)]]], output: OutputDefinition) = accumulatedOutputs flatMap { accumulated =>
      // Extract the valid pairs from the job outputs accumulated so far, and add to it the inputs (outputs can also reference inputs)
      val allKnownValues = accumulated match {
        case Valid(outputs) => outputs.toMap[String, WdlValue] ++ knownValues
        case Invalid(_) => knownValues
      }

      // Attempt to evaluate the expression using all known values
      def evaluateOutputExpression: OutputResult[WdlValue] = {
        EitherT { Try(output.expression.evaluateValue(allKnownValues, ioFunctions)).map(_.toEither) }
      }

      // Attempt to coerce the wdlValue to the desired output type
      def coerceOutputValue(wdlValue: WdlValue, coerceTo: WdlType): OutputResult[WdlValue] = {
        fromEither[Try](
          // TODO WOM: coerceRawValue should return an ErrorOr
          coerceTo.coerceRawValue(wdlValue).toEither.leftMap(t => NonEmptyList.one(t.getMessage))
        )
      }

      /*
        * Go through evaluation, coercion and post processing.
        * Transform the result to a validated Try[ErrorOr[(String, WdlValue)]] with toValidated
        * If we have a valid pair, add it to the previously accumulated outputs, otherwise combine the Nels of errors
       */
      (for {
        evaluated <- evaluateOutputExpression
        coerced <- coerceOutputValue(evaluated, output.womType)
        postProcessed <- fromEither[Try](postMapper(coerced).toErrorOr.toEither)
        pair = output.name -> postProcessed
      } yield pair).toValidated map { evaluatedOutput: ErrorOr[(String, WdlValue)] =>
        (accumulated, evaluatedOutput) mapN {
          case (validValues, validOutput) => validValues :+ validOutput
        }
      }
    }
    
    val emptyValue = Success(List.empty[(String, WdlValue)].validNel): Try[ErrorOr[List[(String, WdlValue)]]]

    // Fold over the outputs to evaluate them in order, map the result to an EvaluatedJobOutputs
    jobDescriptor.call.callable.outputs.foldLeft(emptyValue)(foldFunction) match {
      case Success(Valid(outputs)) => ValidJobOutputs(
        // Wrap the wdlValues in JobOutput objects
        outputs.toMap map {
          case (name, wdlValue) => name -> JobOutput(wdlValue)
        }
      )
      case Success(Invalid(errors)) => InvalidJobOutputs(errors)
      case Failure(exception) => JobOutputsEvaluationException(exception)
    }
  }
}
