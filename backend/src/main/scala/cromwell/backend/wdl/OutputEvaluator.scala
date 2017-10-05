package cromwell.backend.wdl

import cats.Applicative
import cats.data.EitherT._
import cats.data.{EitherT, NonEmptyList}
import cats.instances.try_._
import cats.instances.list._
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
  type OutputResult[A] = EitherT[Try, NonEmptyList[String], A]

  implicit val composedApplicative = Applicative[Try] compose Applicative[ErrorOr]

  def evaluateOutputs(jobDescriptor: BackendJobDescriptor,
                      ioFunctions: IoFunctionSet,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)): OutputResult[CallOutputs] = {
    val knownValues = jobDescriptor.inputDeclarations map {
      case (declaration, value) => declaration.name -> value
    }
    
    def evaluateOutputExpression(expression: WomExpression): OutputResult[WdlValue] = {
      Try(expression.evaluateValue(knownValues, ioFunctions)) match {
        case Success(errorOrValue) => fromEither[Try](errorOrValue.toEither)
        case Failure(ex) => liftT[Try, NonEmptyList[String], WdlValue](Failure(ex))
      }
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
    
    jobDescriptor.call.callable.outputs.traverse[OutputResult, (String, JobOutput)]({ output =>
      validateOutput(output) map { output.name -> _ }
    }).map(_.toMap)
  }
}
