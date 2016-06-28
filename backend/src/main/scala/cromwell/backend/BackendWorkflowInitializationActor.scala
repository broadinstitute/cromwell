package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.backend.wdl.OnlyPureFunctions
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.services.MetadataServiceActor.PutMetadataAction
import cromwell.services.{MetadataEvent, MetadataKey, MetadataValue, ServiceRegistryClient}
import wdl4s.{NoLookup, Task, WdlExpression}
import wdl4s.types._
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

import scala.concurrent.Future
import scala.util.{Failure, Success}

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

}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendWorkflowInitializationActor extends BackendWorkflowLifecycleActor with ServiceRegistryClient with ActorLogging {

  /**
    * Answers the question "does this expression evaluate to a type which matches the predicate".
    *
    * Except the expression may not be evaluable yet, so this is an early sanity check rather than a comprehensive validation.
    */
  protected def wdlTypePredicate(valueRequired: Boolean, predicate: WdlType => Boolean)(wdlExpressionMaybe: Option[WdlExpression]): Boolean = {
    wdlExpressionMaybe match {
      case None => !valueRequired
      case Some(wdlExpression) =>
        wdlExpression.evaluate(NoLookup, OnlyPureFunctions) map (_.wdlType) match {
          case Success(wdlType) => predicate(wdlType)
          case Failure(_) => true // If we can't evaluate it, we'll let it pass for now...
        }
    }
  }

  protected def continueOnReturnCodePredicate(valueRequired: Boolean)(wdlExpressionMaybe: Option[WdlExpression]): Boolean = {
    def isInteger(s: String) = s.forall(_.isDigit) && !s.isEmpty

    wdlExpressionMaybe match {
      case None => !valueRequired
      case Some(wdlExpression) =>
        wdlExpression.evaluate(NoLookup, OnlyPureFunctions) match {
          case Success(wdlValue) => wdlValue match {
            case WdlInteger(_) => true
            case WdlString(_) => isInteger(wdlValue.valueString)
            case WdlArray(WdlArrayType(WdlIntegerType), _) => true
            case WdlArray(WdlArrayType(WdlStringType), elements) => elements forall { x => isInteger(x.valueString) }
            case WdlBoolean(_) => true
            case _ => false
          }
          case Failure(throwable) =>
            throw new RuntimeException(s"Expression evaluation failed due to $throwable: $wdlExpression", throwable)
        }
    }
  }

  protected def runtimeAttributeValidators: Map[String, Option[WdlExpression] => Boolean]

  protected def publishWorkflowRoot(workflowRoot: String) = {
    serviceRegistryActor ! PutMetadataAction(MetadataEvent(MetadataKey(workflowDescriptor.id, None, WorkflowMetadataKeys.WorkflowRoot), MetadataValue(workflowRoot)))
  }

  private def validateRuntimeAttributes: Future[Unit] = {

    def badRuntimeAttrsForTask(task: Task) = {
      runtimeAttributeValidators map { case (attributeName, validator) =>
        val expression = task.runtimeAttributes.attrs.get(attributeName)
        attributeName -> (expression, validator(expression))
      } collect {
        case (name, (expression, false)) => s"Task ${task.name} has an invalid runtime attribute $name = ${expression map {_.valueString} getOrElse "!! NOT FOUND !!"}"
      }
    }

    workflowDescriptor.workflowNamespace.tasks flatMap badRuntimeAttrsForTask match {
      case errors if errors.isEmpty => Future.successful(())
      case errors => Future.failed(new IllegalArgumentException(errors.mkString(". ")))
    }
  }

  def receive: Receive = LoggingReceive {
    case Initialize => performActionThenRespond(initSequence(), onFailure = InitializationFailed)
    case Abort => abortInitialization()
  }

  /**
    * Our predefined sequence to run during preStart
    */
  final def initSequence() = for {
    _ <- validateRuntimeAttributes
    _ <- validate()
    data <- beforeAll()
  } yield InitializationSuccess(data)

  /**
    * Abort all initializations.
    */
  def abortInitialization(): Unit

  /**
    * A call which happens before anything else runs
    */
  def beforeAll(): Future[Option[BackendInitializationData]]

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  def validate(): Future[Unit]

}
