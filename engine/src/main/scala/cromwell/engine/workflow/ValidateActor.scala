package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.backend.{Backend, CromwellBackend}
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.workflow.ValidateActor.{ValidationResult, ValidationFailure, ValidationSuccess}
import cromwell.util.TryUtil
import spray.json
import spray.json.{JsObject, _}
import cromwell.webservice.APIResponse
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.json._
import wdl4s._
import spray.httpx.SprayJsonSupport._
import cromwell.webservice.WorkflowJsonSupport._

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Try, Failure, Success}

object ValidateActor {
  def props(): Props = Props(new ValidateActor())

  sealed trait ValidateActorMessage
  sealed trait ValidationResult extends ValidateActorMessage
  // This ADT contains all the resolved / validated components of the workflow (minus Declarations).
  // It should be used in the future to indicate if the workflow succeeded validation. The callee should reuse these materialized information.
  case class ValidationSuccess(namespaceWithWorkflow: WdlNamespaceWithWorkflow,
                               coercedInputs: Option[WorkflowCoercedInputs],
                               workflowOptions: Option[WorkflowOptions],
                               runtimeAttributes: Seq[Set[String]]) extends ValidationResult
  case class ValidationFailure(reason: Throwable) extends ValidationResult
  case class ValidateWorkflow(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[WdlJson]) extends ValidateActorMessage
}

// TODO: Declarations cannot be validated here currently because of it's dependency on EngineFunctions, and in turn on IOManager and WfContext
class ValidateActor()
  extends Actor with LazyLogging {

  import ValidateActor.ValidateWorkflow
  import context.dispatcher

  override def receive = {
    case ValidateWorkflow(wdlSource, workflowInputs, workflowOptions) =>
      val requester = sender()
        // `validateAll(..)` guarantees we will get a `ValidationResult`, and nothing else
      validateAll(wdlSource, workflowInputs, workflowOptions) map { requester ! _ }
    case unknownMsg => logger.error(s"${this.getClass.getName} received an unknown message: $unknownMsg")
  }

  /**
    * This function collectively validates:
    * 1.) Loading of the namespace,
    * 2.) Raw inputs and it's coercion
    * 3.) Workflow Options
    * 4.) RuntimeAttributes w.r.t. a backend
    * @return The result of the validation as instance of `ValidationResult`
    */
  private def validateAll(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[WdlJson]): Future[ValidationResult] = {
    (for {
      namespaceWithWorkflow <- Future(WdlNamespaceWithWorkflow.load(wdlSource))
      validatedInputs <- Future(validateInputs(workflowInputs, namespaceWithWorkflow))
      validatedWorkflowOptions <- Future(validateWorkflowOptions(workflowOptions, CromwellBackend.backend()))
      validatedRuntimeAttrs <- Future(validateRuntimeAttributes(namespaceWithWorkflow))
    } yield ValidationSuccess(namespaceWithWorkflow, validatedInputs, validatedWorkflowOptions, validatedRuntimeAttrs)) recover {
      case reason => ValidationFailure(reason)
    }
  }

  // TODO: With PBE, this should be defined in the backend.
  private def validateRuntimeAttributes(namespaceWithWorkflow: WdlNamespaceWithWorkflow): Seq[Set[String]] = {
    TryUtil.sequence(namespaceWithWorkflow.workflow.calls map {
      call => CromwellRuntimeAttributes.validateKeys(call.task.runtimeAttributes.attrs.keySet, CromwellBackend.backend().backendType)
    }) match {
      case Success(validatedRuntimeAttrs) => validatedRuntimeAttrs
      case Failure(reason) => throw new IllegalArgumentException("Failed to validate runtime attributes.", reason)
    }
  }

  private def validateInputs(workflowInputs: Option[WdlJson], namespaceWithWorkflow: WdlNamespaceWithWorkflow): Option[WorkflowCoercedInputs] = {
    def rawInputs(json: WdlJson): Map[String, JsValue] = {
      Try(json.parseJson) match {
        case Success(JsObject(inputs)) => inputs
        case Failure(reason: Throwable) => throw new Throwable(reason.getMessage)
        case _ => throw new IllegalArgumentException(s"Failed to parse the workflow inputs JSON: $json")
      }
    }

    workflowInputs match {
      case Some(json) =>
        namespaceWithWorkflow.coerceRawInputs(rawInputs(json)) match {
          case Success(coercedInputs) => Option(coercedInputs)
          case Failure(reason: ThrowableWithErrors) => throw new ValidationException(reason.message, reason.errors)
          case Failure(reason) => throw new IllegalArgumentException("Failed to coerce inputs", reason)
        }
      case _ => None
    }
  }

  // TODO: With PBE, this should be defined in the backend.
  private def validateWorkflowOptions(workflowOptions: Option[WdlJson], backend: Backend): Option[WorkflowOptions] = {
    def validateBackendOptions(options: WorkflowOptions, backend: Backend): WorkflowOptions = {
      try {
        backend.assertWorkflowOptions(options)
        options
      } catch {
        case e: Exception => throw new IllegalArgumentException(s"Workflow has invalid options for backend ${backend.backendType}: ${e.getMessage}", e)
      }
    }

    workflowOptions match {
      case Some(optionsJson) =>
        WorkflowOptions.fromJsonString(optionsJson) match {
          case Success(o) => Option(validateBackendOptions(o, backend))
          case Failure(e) => throw new IllegalArgumentException(s"Failed to validate the workflow options JSON: ${e.getMessage}", e)
        }
      case _ => None
    }
  }
}
