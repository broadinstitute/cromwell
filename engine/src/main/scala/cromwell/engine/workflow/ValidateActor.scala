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
  private val tag = "ValidateActor"

  def props(wdlSource: WdlSource, wdlJson: Option[WdlJson], workflowOptions: Option[WdlJson]): Props = {
    Props(new ValidateActor(wdlSource, wdlJson, workflowOptions))
  }

  sealed trait ValidateActorMessage
  sealed trait ValidationResult extends ValidateActorMessage
  // This ADT contains all the resolved / validated components of the workflow (minus Declarations).
  // It should be used in the future to indicate if the workflow succeeded validation. The callee should reuse these materialized information.
  case class ValidationSuccess(namespaceWithWorkflow: NamespaceWithWorkflow,
                               coercedInputs: Option[WorkflowCoercedInputs],
                               workflowOptions: Option[WorkflowOptions],
                               runtimeAttributes: Seq[Set[String]]) extends ValidationResult
  case class ValidationFailure(reason: Throwable) extends ValidationResult
  case object ValidateWorkflow extends ValidateActorMessage
}

// TODO: Declarations cannot be validated here currently because of it's dependency on EngineFunctions, and in turn on IOManager and WfContext
class ValidateActor(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[String])
  extends Actor with LazyLogging {

  import ValidateActor.{ValidateWorkflow, tag}
  import context.dispatcher

  override def receive = {
    case ValidateWorkflow =>
      validateWorkflow(sender())
    // NOTE: self shuts down when the parent PerRequest shuts down
  }

  private def validateWorkflow(sentBy: ActorRef): Unit = {
    logger.info(s"$tag for $sentBy")
    // TODO: IMO, the responsibility for returning a status code should be a part of the client. This class should only return a validation result,
    // and the callee should convert that to a status code or anything that it wants
    validateAll() map {
      case _: ValidationSuccess =>
        logger.info(s"$tag success $sentBy")
        sentBy ! RequestComplete(
          StatusCodes.OK,
          APIResponse.success("Validation succeeded."))

      case ValidationFailure(reason) =>
        val messageOrBlank = Option(reason.getMessage).mkString
        logger.info(s"$tag error $sentBy: $messageOrBlank")
        sentBy ! RequestComplete(
          StatusCodes.BadRequest,
          APIResponse.fail(reason))
    }
  }

  /**
    * This function collectively validates:
    * 1.) Loading of the namespace,
    * 2.) Raw inputs and it's coercion
    * 3.) Workflow Options
    * 4.) RuntimeAttributes w.r.t. a backend
    * @return
    */
  private def validateAll(): Future[ValidationResult] = {
    (for {
      namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
      validatedInputs <- Future(validateInputs(namespaceWithWorkflow))
      validatedWorkflowOptions <- Future(validateWorkflowOptions(CromwellBackend.backend()))
      validatedRuntimeAttrs <- Future(validateRuntimeAttributes(namespaceWithWorkflow))
    } yield ValidationSuccess(namespaceWithWorkflow, validatedInputs, validatedWorkflowOptions, validatedRuntimeAttrs)) recover {
      case reason => ValidationFailure(reason)
    }
  }

  // TODO: With PBE, this should be defined in the backend.
  private def validateRuntimeAttributes(namespaceWithWorkflow: NamespaceWithWorkflow): Seq[Set[String]] = {
    TryUtil.sequence(namespaceWithWorkflow.workflow.calls map {
      call => CromwellRuntimeAttributes.validateKeys(call.task.runtimeAttributes.attrs.keySet, CromwellBackend.backend().backendType)
    }) match {
      case Success(validatedRuntimeAttrs) => validatedRuntimeAttrs
      case Failure(reason) => throw new IllegalArgumentException("Failed to validate runtime attributes.", reason)
    }
  }

  private def validateInputs(namespaceWithWorkflow: NamespaceWithWorkflow): Option[WorkflowCoercedInputs] = {
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
  private def validateWorkflowOptions(backend: Backend): Option[WorkflowOptions] = {
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
