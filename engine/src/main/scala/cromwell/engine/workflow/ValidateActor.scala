package cromwell.engine.workflow

import akka.actor.{Actor, Props}
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.{Backend, CromwellBackend}
import cromwell.engine.ErrorOr
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.backend.{Backend, CromwellBackend}
import cromwell.engine.workflow.ValidateActor.{ValidationFailure, ValidationResult, ValidationSuccess}
import cromwell.util.TryUtil
import spray.json._
import wdl4s._

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._

object ValidateActor {
  private val EmptyJson: String = "{}"

  def props(): Props = Props(new ValidateActor())

  sealed trait ValidateActorMessage
  sealed trait ValidationResult extends ValidateActorMessage

  // This ADT contains all the resolved / validated components of the workflow (minus Declarations).
  // It should be used in the future to indicate if the workflow succeeded validation. The callee should reuse these materialized information.
  case class ValidationSuccess(namespaceWithWorkflow: NamespaceWithWorkflow,
                               coercedInputs: Option[WorkflowCoercedInputs],
                               workflowOptions: Option[WorkflowOptions],
                               runtimeAttributes: Seq[Set[String]]) extends ValidationResult
  case class ValidationFailure(reason: Throwable) extends ValidationResult
  case class ValidateWorkflow(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[WdlJson]) extends ValidateActorMessage
}

// TODO: Declarations cannot be validated here currently because of it's dependency on EngineFunctions, and in turn on IOManager and WfContext
class ValidateActor() extends Actor with LazyLogging {

  import ValidateActor.{EmptyJson, ValidateWorkflow}
  import context.dispatcher

  override def receive = {
    case ValidateWorkflow(wdlSource, workflowInputs, workflowOptions) =>
      val requester = sender()
      val backend = workflowOptions map { CromwellBackend.getBackendFromOptions } getOrElse CromwellBackend.defaultBackend
      // `validateAll(..)` guarantees we will get a `ValidationResult`, and nothing else
      validateAll(wdlSource, workflowInputs, workflowOptions, backend) map {
        requester ! _
      }
    case unknownMsg => logger.error(s"${this.getClass.getName} received an unknown message: $unknownMsg")
  }

  /**
    * This function collectively validates:
    * 1.) Loading of the namespace
    * 2.) Raw inputs and their coercions
    * 3.) Workflow Options
    * 4.) RuntimeAttributes w.r.t. a backend
    *
    * @return The result of the validation as instance of `ValidationResult`
    */
  private def validateAll(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[WdlJson], backend: Backend): Future[ValidationResult] = Future {
    val namespace = validateNamespace(wdlSource)
    val options = validateWorkflowOptions(backend, workflowOptions: Option[WdlJson])
    val result = (namespace |@| options) { (_, _) } flatMap {
      case (ns, opt) =>
        val inputs = validateInputs(ns, workflowInputs)
        val runtimeAttributes = validateRuntimeAttributes(ns, backend)
        (inputs |@| runtimeAttributes) ((in, runAttr) => ValidationSuccess(ns, Option(in), Option(opt), runAttr))
    }

    result match {
      case scalaz.Success(success) => success
      case scalaz.Failure(failures) =>
        ValidationFailure(new ValidationException(s"Workflow input processing failed.", failures))
    }
  }

  private def validateNamespace(source: WdlSource): ErrorOr[NamespaceWithWorkflow] = {
    try {
      NamespaceWithWorkflow.load(source).successNel
    } catch {
      case e: Exception => s"Unable to load namespace from workflow: ${e.getMessage}".failureNel
    }
  }

  private def validateInputs(namespaceWithWorkflow: NamespaceWithWorkflow, workflowInputs: Option[WdlJson]): ErrorOr[WorkflowCoercedInputs] = {
    def validateRawInputs(json: WdlJson): ErrorOr[Map[String, JsValue]] = {
      Try(json.parseJson) match {
        case Success(JsObject(inputs)) => inputs.successNel
        case Failure(reason: Throwable) => s"Workflow contains invalid inputs JSON: ${reason.getMessage}".failureNel
        case _ => s"Workflow inputs JSON cannot be parsed to JsObject: $json".failureNel
      }
    }

    validateRawInputs(workflowInputs.getOrElse(EmptyJson)) match {
      case scalaz.Success(inputs) =>
        namespaceWithWorkflow.coerceRawInputs(inputs) match {
          case Success(r) => r.successNel
          case Failure(e: ThrowableWithErrors) => scalaz.Failure(e.errors)
          case Failure(e) => e.getMessage.failureNel
        }
      case scalaz.Failure(f) => f.failure
    }
  }

  // TODO: With PBE, this should be defined in the backend.
  private def validateWorkflowOptions(backend: Backend, workflowOptions: Option[WdlJson]): ErrorOr[WorkflowOptions] = {
    def validateBackendOptions(backend: Backend, workflowOpt: WorkflowOptions): ErrorOr[WorkflowOptions] = {
      try {
        backend.assertWorkflowOptions(workflowOpt)
        workflowOpt.successNel
      } catch {
        case e: Exception => s"Workflow has invalid options for backend ${backend.backendType}: ${e.getMessage}".failureNel
      }
    }

    WorkflowOptions.fromJsonString(workflowOptions.getOrElse(EmptyJson)) match {
      case Success(o) => validateBackendOptions(backend, o)
      case Failure(e) => s"Workflow contains invalid options JSON: ${e.getMessage}".failureNel
    }
  }

  // TODO: With PBE, this should be defined in the backend.
  // TODO: Add CromwellRuntimeAttributes as a dependency for this actor (arg in ctor) in case is not moved to specific
  // backend when PBE is merged.
  private def validateRuntimeAttributes(namespaceWithWorkflow: NamespaceWithWorkflow, backend: Backend): ErrorOr[Seq[Set[String]]] = {
    TryUtil.sequence(namespaceWithWorkflow.workflow.calls map {
      call => CromwellRuntimeAttributes.validateKeys(call.task.runtimeAttributes.attrs.keySet, backend.backendType)
    }) match {
      case Success(validatedRuntimeAttrs) => validatedRuntimeAttrs.successNel
      case Failure(reason) => "Failed to validate runtime attributes.".failureNel
    }
  }
}
