package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.backend.validation.TryUtils
import cromwell.backend.validation.exception.ValidationAggregatedException
import cromwell.core._
import lenthall.exception.AggregatedException
import wdl4s.expression.{NoFunctions, ValueEvaluator}
import wdl4s.{Call, WdlExpression}

import scala.concurrent.Future
import scalaz.{Validation, Failure, NonEmptyList}

object BackendWorkflowInitializationActor {

  // Commands
  sealed trait BackendWorkflowInitializationActorCommand extends BackendWorkflowLifecycleActorCommand

  case object Initialize extends BackendWorkflowInitializationActorCommand

  final case class Abort(jobKey: BackendJobDescriptorKey) extends BackendWorkflowInitializationActorCommand

  // Responses
  sealed trait BackendWorkflowInitializationActorResponse extends BackendWorkflowLifecycleActorResponse

  sealed trait InitializationResponse extends BackendWorkflowInitializationActorResponse

  case object InitializationSuccess extends InitializationResponse

  case class InitializationFailed(reason: Throwable) extends Exception with InitializationResponse

}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendWorkflowInitializationActor extends BackendWorkflowLifecycleActor with ActorLogging {

  def receive: Receive = LoggingReceive {
    case Initialize => performActionThenRespond(initSequence(), onFailure = InitializationFailed)
    case AbortWorkflow => performActionThenRespond(abortInitialization(), onFailure = BackendWorkflowAbortFailedResponse)
  }

  /**
    * Our predefined sequence to run during preStart
    */
  final def initSequence() = for {
    _ <- validate()
    _ <- beforeAll()
  } yield InitializationSuccess

  /**
    * Abort all initializations.
    */
  def abortInitialization(): Future[WorkflowAbortResponse]

  /**
    * A call which happens before anything else runs
    */
  def beforeAll(): Future[Unit]

  /**
    * Validates runtime attributes for one specific call.
    *
    * @param runtimeAttributes Runtime Attributes with already evaluated values.
    * @return If all entries from runtime attributes section are valid Success otherwise
    *         Failure with the aggregation of errors.
    */
  def validateRuntimeAttributes(runtimeAttributes: EvaluatedRuntimeAttributes): Future[scalaz.Validation[NonEmptyList[String], Unit]]

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  protected def validate(): Future[Unit] = {
    def eval(call: Call) = validateRuntimeAttributes(evaluateRuntimeAttributesWdlExpressions(call))

    val validationResult: Future[Seq[(Call, Validation[NonEmptyList[String], Unit])]] =
      Future.sequence(calls map { call => eval(call) map { valRes => (call, valRes) } })

    val failedValidations = validationResult.map(evaluateValidationResultPerCall)

    //TODO: Add recover or not but refactor this as part of https://github.com/broadinstitute/cromwell/issues/725.
    failedValidations map { failures =>
      if (failures.nonEmpty) {
        throw new AggregatedException("Runtime attribute validation failed", failures)
      }
    }
  }

  /**
    * Creates WDL expression evaluator for the specific call.
    *
    * @param call Call which contains WDL expressions in it.
    * @return A WDL ValueEvaluator.
    */
  protected def createEvaluator(call: Call): ValueEvaluator = {
    val declarations = workflowDescriptor.workflowNamespace.workflow.declarations ++ call.task.declarations
    val knownInputs = workflowDescriptor.inputs
    val lookup = WdlExpression.standardLookupFunction(knownInputs, declarations, NoFunctions)
    new ValueEvaluator(lookup, NoFunctions)
  }

  /**
    * Evaluates WDL expressions within call runtime attributes.
    *
    * @param call Call which contains WDL expressions in it.
    * @return Evaluated runtime attributes.
    */
  private def evaluateRuntimeAttributesWdlExpressions(call: Call): EvaluatedRuntimeAttributes = {
    val evaluator = createEvaluator(call)
    def evaluate(wdlExpression: WdlExpression) = evaluator.evaluate(wdlExpression.ast)
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues evaluate
    TryUtils.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
  }

  /**
    * Checks if there was any validation failure.
    *
    * @param evaluation Tuple with a call and related validation result.
    * @return A list with all validation exceptions.
    */
  private def evaluateValidationResultPerCall(evaluation: Seq[(Call, ErrorOr[Unit])]): Seq[ValidationAggregatedException] = {
    evaluation collect {
      case (call, validation: Failure[NonEmptyList[String]]) =>
        validation match {
          case failureWithStringErr: Failure[NonEmptyList[String]] =>
            ValidationAggregatedException(s"Runtime attribute validation failed for task '${call.taskFqn}'.", failureWithStringErr.e.list)
        }
    }
  }
}
