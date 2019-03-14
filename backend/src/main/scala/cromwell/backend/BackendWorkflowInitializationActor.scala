package cromwell.backend

import akka.actor.{ActorLogging, ActorRef}
import akka.event.LoggingReceive
import cats.data.Validated.{Invalid, Valid, _}
import cats.data.{NonEmptyList, ValidatedNel}
import cats.instances.list._
import cats.syntax.traverse._
import common.collections.EnhancedCollections._
import common.validation.Validation._
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.backend.async.{RuntimeAttributeValidationFailure, RuntimeAttributeValidationFailures}
import cromwell.backend.validation.ContinueOnReturnCodeValidation
import cromwell.core._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import mouse.all._
import wom.expression.{NoIoFunctionSet, WomExpression}
import wom.graph.CommandCallNode
import wom.types.WomType
import wom.values.WomValue

import scala.concurrent.Future
import scala.util.Try

object BackendWorkflowInitializationActor {

  // Commands
  sealed trait BackendWorkflowInitializationActorCommand extends BackendWorkflowLifecycleActorCommand
  case object Initialize extends BackendWorkflowInitializationActorCommand
  case object Abort extends BackendWorkflowInitializationActorCommand

  // Responses
  sealed trait BackendWorkflowInitializationActorResponse extends BackendWorkflowLifecycleActorResponse
  sealed trait InitializationResponse extends BackendWorkflowInitializationActorResponse
  case class InitializationSuccess(backendInitializationData: Option[BackendInitializationData]) extends InitializationResponse
  case class InitializationFailed(reason: Throwable) extends Exception with InitializationResponse

  /**
    * This method calls into `runtimeAttributeValidators` to validate runtime attribute expressions.
    * The initialization-time validators try to fully `evaluate` (not `evaluateType`) the runtime attribute expression.
    * Because these validators do full evaluations but at initialization time we don't have the actual call inputs, it's
    * possible a legitimate runtime attribute expression evaluation could fail to evaluate.  These validators must
    * therefore give runtime attribute expressions which failed to evaluate the benefit of the doubt and pass them.
    *
    * Perhaps a better way of doing initialization-time runtime attribute validation:
    *
    * - Currently `default_runtime_attributes` are values and not expressions.  These should all be typechecked
    *   irrespective of whether there is actual fallback to using their values due to missing task attributes.
    *
    * - It appears it should be possible to fully typecheck runtime attribute expressions at workflow initialization
    *   time if task (not call!) declarations and inputs are considered.  We have all tasks at workflow initialization
    *   time, and the task declarations know their types.  We won't actually have a call until all of its inputs are
    *   available, but tasks are available in the `WdlNamespace`.  A task-aware type evaluation might look like:
    *
    *   {{{
    *     val task: Task = ...
    *     val expression: WdlExpression = ...
    *     def lookup(name: String): WomType = task.declarations.find(_.name == name).get.womType
    *     val expressionType = expression.evaluateType(lookup, OnlyPureFunctions)
    *   }}}
    *
    * - Assuming the previous two checks pass, an attempt should be made to evaluate the expression.  It's possible
    *   the expression could be a function of call inputs or declarations which are unavailable at initialization time,
    *   so these failed evaluations must continue to be accepted as possibly correct expressions.
    *
    * - Assuming the expression evaluation succeeds, it should be possible to do real business logic validation.
    *
    * - It would be nice to memoize as much of the work that gets done here as possible so it doesn't have to all be
    *   repeated when the various `FooRuntimeAttributes` classes are created, in the spirit of #1076.
    */
    def validateRuntimeAttributes(
                                   taskName: String,
                                   defaultRuntimeAttributes: Map[String, WomValue],
                                   runtimeAttributes: Map[String, WomExpression],
                                   runtimeAttributeValidators: Map[String, Option[WomExpression] => Boolean]
                                 ): ValidatedNel[RuntimeAttributeValidationFailure, Unit] = {

      //This map append will overwrite default key/values with runtime settings upon key collisions
      val lookups = defaultRuntimeAttributes.safeMapValues(_.asWomExpression) ++ runtimeAttributes

      runtimeAttributeValidators.toList.traverse{
        case (attributeName, validator) =>
          val runtimeAttributeValue: Option[WomExpression] = lookups.get(attributeName)
          validator(runtimeAttributeValue).fold(
            validNel(()),
            Invalid(NonEmptyList.of(RuntimeAttributeValidationFailure(taskName, attributeName, runtimeAttributeValue)))
          )
      }.map(_ => ())
  }
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendWorkflowInitializationActor extends BackendWorkflowLifecycleActor with ActorLogging {
  def serviceRegistryActor: ActorRef

  def calls: Set[CommandCallNode]

  /**
    * This method is meant only as a "pre-flight check" validation of runtime attribute expressions during workflow
    * initialization.  This attempts to evaluate an expression without lookups and with only pure functions.  Successful
    * evaluations will be typechecked, failed expression evaluations are assumed (perhaps wrongly) to correspond to
    * valid expressions.  This assumption is made since at the time this method executes it does not have actual call
    * inputs which might be required for evaluation of a runtime attribute expression.
    *
    * Currently this does not consider task declarations or inputs.  So runtime attributes which are functions of task
    * declarations will fail evaluation and return `true` from this predicate, even if the type could be determined
    * to be wrong with consideration of task declarations or inputs.
    */
  protected def womTypePredicate(valueRequired: Boolean, predicate: WomType => Boolean)(womExpressionMaybe: Option[WomExpression]): Boolean = {
    womExpressionMaybe match {
      case None => !valueRequired
      case Some(womExpression: WomExpression) =>
        womExpression.evaluateValue(Map.empty, NoIoFunctionSet) map (_.womType) match {
          case Valid(womType) => predicate(womType)
          case Invalid(_) => true // If we can't evaluate it, we'll let it pass for now...
        }
    }
  }

  /**
    * This predicate is only appropriate for validation during workflow initialization.  The logic does not differentiate
    * between evaluation failures due to missing call inputs or evaluation failures due to malformed expressions, and will
    * return `true` in both cases.
    */
  protected def continueOnReturnCodePredicate(valueRequired: Boolean)(womExpressionMaybe: Option[WomValue]): Boolean = {
    ContinueOnReturnCodeValidation.default(configurationDescriptor.backendRuntimeAttributesConfig).validateOptionalWomValue(womExpressionMaybe)
  }

  protected def runtimeAttributeValidators: Map[String, Option[WomExpression] => Boolean]

  // FIXME: If a workflow executes jobs using multiple backends,
  // each backend will try to write its own workflow root and override any previous one.
  // They should be structured differently or at least be prefixed by the backend name
  protected def publishWorkflowRoot(workflowRoot: String): Unit = {
    serviceRegistryActor ! PutMetadataAction(MetadataEvent(MetadataKey(workflowDescriptor.id, None, WorkflowMetadataKeys.WorkflowRoot), MetadataValue(workflowRoot)))
  }

  protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WomValue]]


  def receive: Receive = LoggingReceive {
    case Initialize => performActionThenRespond(initSequence(), onFailure = InitializationFailed)
    case Abort => abortInitialization()
  }

  /**
    * Our predefined sequence to run during preStart
    */
  final def initSequence(): Future[InitializationSuccess] =
    for {
      defaultRuntimeAttributes <- coerceDefaultRuntimeAttributes(workflowDescriptor.workflowOptions) |> Future.fromTry _
      taskList = calls.toList.map(_.callable).map(t => t.name -> t.runtimeAttributes.attributes)
      _ <- taskList.
            traverse{
              case (name, runtimeAttributes) => validateRuntimeAttributes(name, defaultRuntimeAttributes, runtimeAttributes, runtimeAttributeValidators)
            }.toFuture(errors => RuntimeAttributeValidationFailures(errors.toList))
      _ <- validate()
      data <- beforeAll()
    } yield InitializationSuccess(data)

  /**
    * Abort all initializations.
    */
  def abortInitialization(): Unit = {}

  /**
    * A call which happens before anything else runs
    */
  def beforeAll(): Future[Option[BackendInitializationData]]

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  def validate(): Future[Unit]

}
