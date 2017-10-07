package cromwell.backend.wdl

import cats.Applicative
import cats.data.EitherT._
import cats.data.Validated.{Invalid, Valid}
import cats.data.{EitherT, NonEmptyList}
import cats.instances.list._
import cats.instances.try_._
import cats.syntax.either._
import cats.syntax.traverse._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.{CallOutputs, JobOutput}
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import wdl.types.WdlType
import wdl.values.WdlValue
import wom.callable.Callable.OutputDefinition
import wom.expression.{IoFunctionSet, WomExpression}

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
    val knownValues = jobDescriptor.inputDeclarations map {
      case (declaration, value) => declaration.name -> value
    }
    
    def evaluateOutputExpression(expression: WomExpression): OutputResult[WdlValue] = {
      EitherT { Try(expression.evaluateValue(knownValues, ioFunctions)).map(_.toEither) }
    }
    
    def coerceOutputValue(wdlValue: WdlValue, coerceTo: WdlType): OutputResult[WdlValue] = {
      fromEither[Try](
        coerceTo.coerceRawValue(wdlValue).toEither
        .leftMap(t => NonEmptyList.one(t.getMessage))
      )
    }
    
    def validateOutput(output: OutputDefinition): OutputResult[JobOutput] = for {
      evaluated <- evaluateOutputExpression(output.expression)
      coerced <- coerceOutputValue(evaluated, output.womType)
      postProcessed <- fromEither[Try](postMapper(coerced).toErrorOr.toEither)
    } yield JobOutput(postProcessed)
    
    type TryErrorOr[A] = Try[ErrorOr[A]]
    
    /*
      * Traverse the list of outputs and evaluate, coerce and post process each of them.
     */
    jobDescriptor.call.callable.outputs.traverse[TryErrorOr, (String, JobOutput)]({ output =>
      // Validate the output and get an OutputResult
      validateOutput(output)
        // Create a tuple of output name -> value
        .map(output.name -> _)
        // Access the value of the OutputResult which is a Try[Either[NonEmptyList[String], (String, JobOutput)]]
        .value
        // map the Either in the Try to an ErrorOr to combine all errors using the applicative property of ErrorOr
        .map(_.toValidated)
    })(Applicative[Try] compose Applicative[ErrorOr]) match {
      case Success(Valid(outputs)) => ValidJobOutputs(outputs.toMap)
      case Success(Invalid(errors)) => InvalidJobOutputs(errors)
      case Failure(exception) => JobOutputsEvaluationException(exception)
    }
  }
}
