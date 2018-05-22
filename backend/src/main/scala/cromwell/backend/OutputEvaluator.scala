package cromwell.backend.wdl

import cats.data.EitherT._
import cats.data.Validated.{Invalid, Valid}
import cats.data.{EitherT, NonEmptyList}
import cats.instances.try_._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.util.TryUtil
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CallOutputs
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.{ExpressionBasedOutputPort, OutputPort}
import wom.types.WomType
import wom.values.WomValue

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object OutputEvaluator {
  sealed trait EvaluatedJobOutputs
  case class ValidJobOutputs(outputs: CallOutputs) extends EvaluatedJobOutputs
  case class InvalidJobOutputs(errors: NonEmptyList[String]) extends EvaluatedJobOutputs
  case class JobOutputsEvaluationException(exception: Throwable) extends EvaluatedJobOutputs

  type OutputResult[A] = EitherT[Try, NonEmptyList[String], A]

  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      ioFunctions: IoFunctionSet,
                      postMapper: WomValue => Try[WomValue] = v => Success(v))(implicit ec: ExecutionContext): Future[EvaluatedJobOutputs] = {
    val taskInputValues: Map[String, WomValue] = jobDescriptor.localInputs

    def foldFunction(accumulatedOutputs: Try[ErrorOr[List[(OutputPort, WomValue)]]], output: ExpressionBasedOutputPort) = accumulatedOutputs flatMap { accumulated =>
      // Extract the valid pairs from the job outputs accumulated so far, and add to it the inputs (outputs can also reference inputs)
      val allKnownValues: Map[String, WomValue] = accumulated match {
        case Valid(outputs) =>
          // The evaluateValue methods needs a Map[String, WomValue], use the output port name for already computed outputs
          outputs.toMap[OutputPort, WomValue].map({ case (port, value) => port.internalName -> value }) ++ taskInputValues
        case Invalid(_) => taskInputValues
      }

      // Attempt to evaluate the expression using all known values
      def evaluateOutputExpression: OutputResult[WomValue] =
        EitherT.fromEither[Try] {
          output.expression.evaluateValue(allKnownValues, ioFunctions).toEither
        }

      // Attempt to coerce the womValue to the desired output type
      def coerceOutputValue(womValue: WomValue, coerceTo: WomType): OutputResult[WomValue] = {
        fromEither[Try](
          // TODO WOM: coerceRawValue should return an ErrorOr
          coerceTo.coerceRawValue(womValue).toEither.leftMap(t => NonEmptyList.one(t.getClass.getSimpleName + ": " + t.getMessage))
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

      def enhanceErrorText(s: String): String = s"Bad output '${output.name}': $s"
      evaluated.leftMap(_ map enhanceErrorText).toValidated map { evaluatedOutput: ErrorOr[(OutputPort, WomValue)] =>
        (accumulated, evaluatedOutput) mapN { _ :+ _ }
      }
    }

    val emptyValue = Success(List.empty[(OutputPort, WomValue)].validNel): Try[ErrorOr[List[(OutputPort, WomValue)]]]

    // Fold over the outputs to evaluate them in order, map the result to an EvaluatedJobOutputs
    def fromOutputPorts: EvaluatedJobOutputs = jobDescriptor.taskCall.expressionBasedOutputPorts.foldLeft(emptyValue)(foldFunction) match {
      case Success(Valid(outputs)) => ValidJobOutputs(CallOutputs(outputs.toMap))
      case Success(Invalid(errors)) => InvalidJobOutputs(errors)
      case Failure(exception) => JobOutputsEvaluationException(exception)
    }
    
    /*
      * Because Cromwell doesn't trust anyone, if custom evaluation is provided,
      * still make sure that all the output ports have been filled with values
     */
    def validateCustomEvaluation(outputs: Map[OutputPort, WomValue]): EvaluatedJobOutputs = {
      def toError(outputPort: OutputPort) = s"Missing output value for ${outputPort.identifier.fullyQualifiedName.value}"

      jobDescriptor.taskCall.expressionBasedOutputPorts.diff(outputs.keySet.toList) match {
        case Nil =>
          val errorMessagePrefix = "Error applying postMapper in short-circuit output evaluation"
          TryUtil.sequenceMap(outputs map { case (k, v) => (k, postMapper(v))}, errorMessagePrefix) match {
            case Failure(e) => InvalidJobOutputs(NonEmptyList.of(e.getMessage, e.getStackTrace.take(5).map(_.toString):_*))
            case Success(postMappedOutputs) => ValidJobOutputs(CallOutputs(postMappedOutputs))
          }
        case head :: tail => InvalidJobOutputs(NonEmptyList.of(toError(head), tail.map(toError): _*))
      }
    }

    /*
      * See if the task definition has "short-circuit" for the default output evaluation.
      * In the case of CWL for example, this gives a chance to look for cwl.output.json and use it as the output of the tool,
      * instead of the default behavior of going over each output port of the task and evaluates their expression.
      * If the "customOutputEvaluation" returns None (which will happen if the cwl.output.json is not there, as well as for all WDL workflows),
      * we fallback to the default behavior.
     */
    jobDescriptor.taskCall.customOutputEvaluation(taskInputValues, ioFunctions, ec).value
      .map({
        case Some(Right(outputs)) => validateCustomEvaluation(outputs)
        case Some(Left(errors)) => InvalidJobOutputs(errors)
        // If it returns an empty value, fallback to canonical output evaluation
        case None => fromOutputPorts
      })
  }
}
