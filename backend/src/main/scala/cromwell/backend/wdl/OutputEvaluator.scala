package cromwell.backend.wdl

import cats.data.EitherT._
import cats.data.Validated.{Invalid, Valid}
import cats.data.{EitherT, NonEmptyList}
import cats.instances.try_._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import cromwell.backend.BackendJobDescriptor
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.CallOutputs
import cromwell.backend.io.GlobFunctions
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.{ExpressionBasedOutputPort, OutputPort}
import wom.types.WomType
import wom.values.{WomSingleFile, WomValue}

import scala.util.{Failure, Success, Try}
object OutputEvaluator {
  sealed trait EvaluatedJobOutputs
  case class ValidJobOutputs(outputs: CallOutputs) extends EvaluatedJobOutputs
  case class InvalidJobOutputs(errors: NonEmptyList[String]) extends EvaluatedJobOutputs
  case class JobOutputsEvaluationException(exception: Throwable) extends EvaluatedJobOutputs

  type OutputResult[A] = EitherT[Try, NonEmptyList[String], A]

  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      ioFunctions: IoFunctionSet,
                      postMapper: WomValue => Try[WomValue] = v => Success(v)): EvaluatedJobOutputs = {
    val knownValues: Map[String, WomValue] = jobDescriptor.localInputs
    
    def foldFunction(accumulatedOutputs: Try[ErrorOr[List[(OutputPort, WomValue)]]], output: ExpressionBasedOutputPort) = accumulatedOutputs flatMap { accumulated =>
      // Extract the valid pairs from the job outputs accumulated so far, and add to it the inputs (outputs can also reference inputs)
      val allKnownValues: Map[String, WomValue] = accumulated match {
        case Valid(outputs) => 
          // The evaluateValue methods needs a Map[String, WomValue], use the output port name for already computed outputs
          outputs.toMap[OutputPort, WomValue].map({ case (port, value) => port.name -> value }) ++ knownValues
        case Invalid(_) => knownValues
      }

      // Attempt to evaluate the expression using all known values
      def evaluateOutputExpression: OutputResult[WomValue] = EitherT{ Try{
        output.expression.evaluateValue(allKnownValues, ioFunctions).toEither.map{

          //If the string contained in a filename is actually a glob,
          // The output of this expression is expected to be glob-md5(fileName).list
          case WomSingleFile(fileName) if fileName.contains("*") =>
            WomSingleFile(GlobFunctions.globName(fileName) + ".list")

          case other => other
        }
      }}

      // Attempt to coerce the womValue to the desired output type
      def coerceOutputValue(womValue: WomValue, coerceTo: WomType): OutputResult[WomValue] = {
        fromEither[Try](
          // TODO WOM: coerceRawValue should return an ErrorOr
          coerceTo.coerceRawValue(womValue).toEither.leftMap(t => NonEmptyList.one(t.getMessage))
        )
      }

      /*
        * Go through evaluation, coercion and post processing.
        * Transform the result to a validated Try[ErrorOr[(String, WomValue)]] with toValidated
        * If we have a valid pair, add it to the previously accumulated outputs, otherwise combine the Nels of errors
       */
      val evaluated = for {
        evaluated <- evaluateOutputExpression
        coerced <- coerceOutputValue(evaluated, output.womType)
        postProcessed <- EitherT { postMapper(coerced).map(_.validNelCheck) }: OutputResult[WomValue]
        pair = output -> postProcessed
      } yield pair

      evaluated.toValidated map { evaluatedOutput: ErrorOr[(OutputPort, WomValue)] =>
        (accumulated, evaluatedOutput) mapN { _ :+ _ }
      }
    }
    
    val emptyValue = Success(List.empty[(OutputPort, WomValue)].validNel): Try[ErrorOr[List[(OutputPort, WomValue)]]]

    // Fold over the outputs to evaluate them in order, map the result to an EvaluatedJobOutputs
    jobDescriptor.call.expressionBasedOutputPorts.foldLeft(emptyValue)(foldFunction) match {
      case Success(Valid(outputs)) => ValidJobOutputs(CallOutputs(outputs.toMap))
      case Success(Invalid(errors)) => InvalidJobOutputs(errors)
      case Failure(exception) => JobOutputsEvaluationException(exception)
    }
  }
}
