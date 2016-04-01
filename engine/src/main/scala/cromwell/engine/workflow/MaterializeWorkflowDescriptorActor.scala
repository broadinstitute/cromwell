package cromwell.engine.workflow

import java.nio.file.Paths

import akka.actor.{Actor, Props}
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.{OptionNotFoundException, ErrorOr, WorkflowOptions, WorkflowId}
import cromwell.engine._
import cromwell.engine.analysis.BackendSelector
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.backend._
import cromwell.util.TryUtil
import spray.json.{JsObject,_}
import wdl4s._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._
import lenthall.config.ScalaConfig.EnhancedScalaConfig


object MaterializeWorkflowDescriptorActor {
  def props(): Props = Props(new MaterializeWorkflowDescriptorActor)

  sealed trait MaterializeWorkflowDescriptorActorMessage
  case class MaterializeWorkflow(id: WorkflowId,
                                 workflowSourceFiles: WorkflowSourceFiles,
                                 conf: Config = ConfigFactory.load) extends MaterializeWorkflowDescriptorActorMessage
  sealed trait MaterializationResult extends MaterializeWorkflowDescriptorActorMessage
  case class MaterializationSuccess(workflowDescriptor: WorkflowDescriptor) extends MaterializationResult
  case class MaterializationFailure(reason: Throwable) extends MaterializationResult

  import lenthall.config.ScalaConfig._

  private val DefaultCallCachingValue = false
  private val DefaultLookupDockerHash = false
  private val DefaultWorkflowFailureMode = NoNewCalls.toString

  def configCallCaching(conf: Config) = lookupBooleanWithDefault(conf, "call-caching", "enabled", DefaultCallCachingValue)
  def lookupDockerHash(conf: Config) = lookupBooleanWithDefault(conf, "call-caching", "lookup-docker-hash", DefaultLookupDockerHash)

  private def lookupBooleanWithDefault(conf: Config, stanza: String, key: String, default: Boolean) = {
    (for {
      config <- conf.getConfigOption(stanza)
      value <- config.getBooleanOption(key)
    } yield value) getOrElse default
  }

  def workflowLogOptions(conf: Config): Option[WorkflowLogOptions] = {
    for {
      workflowConfig <- conf.getConfigOption("workflow-options")
      dir <- workflowConfig.getStringOption("workflow-log-dir") if !dir.isEmpty
      temporary <- workflowConfig.getBooleanOption("workflow-log-temporary") orElse Option(true)
    } yield WorkflowLogOptions(Paths.get(dir), temporary)
  }
}

class MaterializeWorkflowDescriptorActor() extends Actor with LazyLogging {

  import MaterializeWorkflowDescriptorActor._

  override def receive = {
    case MaterializeWorkflow(workflowId, workflowSourceFiles, conf) =>
      buildWorkflowDescriptor(workflowId, workflowSourceFiles, conf) match {
        case scalaz.Success(descriptor) => sender() ! MaterializationSuccess(descriptor)
        case scalaz.Failure(error) => sender() ! MaterializationFailure(
          new IllegalArgumentException() with ThrowableWithErrors {
            val message = s"Workflow input processing failed."
            val errors = error
          })
      }
    case unknownMsg => logger.error(s"${this.getClass.getName} received an unknown message: $unknownMsg")
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFiles,
                                      conf: Config): ErrorOr[WorkflowDescriptor] = {

    def buildWorkflowDescriptorInner(namespace: NamespaceWithWorkflow,
                                             workflowOptions: WorkflowOptions,
                                             rawInputs: Map[String, JsValue],
                                             defaultBackend: Backend,
                                             backendAssignments: Map[String, Backend],
                                             workflowFailureMode: WorkflowFailureMode): ErrorOr[WorkflowDescriptor] = {
      validateCoercedInputs(rawInputs, namespace) flatMap { coercedInputs =>
        val workflowRootPath = defaultBackend.buildWorkflowRootPath(defaultBackend.rootPath(workflowOptions), namespace.workflow.unqualifiedName, id)
        val wfContext = new WorkflowContext(workflowRootPath)
        val fileSystems = defaultBackend.fileSystems(workflowOptions)
        val engineFunctions = defaultBackend.engineFunctions(fileSystems, wfContext)

        validateDeclarations(namespace, workflowOptions, coercedInputs, engineFunctions) flatMap { declarations =>
          WorkflowDescriptor(id, sourceFiles, workflowOptions, workflowLogOptions(conf), rawInputs, namespace, coercedInputs, declarations, defaultBackend,
            backendAssignments, configCallCaching(conf), lookupDockerHash(conf), workflowFailureMode, wfContext, fileSystems).successNel
        }
      }
    }

    val namespaceValidation = validateNamespace(sourceFiles.wdlSource)
    val workflowOptionsValidation = validateWorkflowOptions(sourceFiles.workflowOptionsJson)
    (namespaceValidation |@| workflowOptionsValidation) {
      (_, _)
    } flatMap { case (namespace, workflowOptions) =>
      val defaultBackend = CromwellBackend.getDefaultBackendForWorkflowAfterConsideringWorkflowOptions(workflowOptions)
      val backendAssignmentsValidation = materializeValidManifestedBackendAssignments(namespace, defaultBackend)
      val rawInputsValidation = validateRawInputs(sourceFiles.inputsJson)
      val failureModeValidation = validateWorkflowFailureMode(workflowOptions, conf)
      (backendAssignmentsValidation |@| rawInputsValidation |@| failureModeValidation) {
        (_, _, _)
      } flatMap { case (backendAssignments, rawInputs, failureMode) =>
        validateRuntimeAttributes(backendAssignments, namespace) flatMap { _ =>
          buildWorkflowDescriptorInner(namespace, workflowOptions, rawInputs, defaultBackend, backendAssignments, failureMode)
        }
      }
    }
  }

  private def materializeValidManifestedBackendAssignments(namespace: NamespaceWithWorkflow, defaultBackend: Backend): ErrorOr[Map[FullyQualifiedName, Backend]] = {
    // TODO: A Try that's so ugly even its parents would struggle to love it. Crying out for better validation:
    Try(namespace.workflow.calls map { call =>
      (call.fullyQualifiedName, BackendSelector.selectBackend(Option(defaultBackend), call).get)
    } toMap) match {
      case Success(assignments) => assignments.successNel
      case Failure(t) => t.getMessage.failureNel
    }
  }

  private def validateDeclarations(namespace: NamespaceWithWorkflow,
                                   options: WorkflowOptions,
                                   coercedInputs: WorkflowCoercedInputs,
                                   engineFunctions: WorkflowEngineFunctions): ErrorOr[WorkflowCoercedInputs] = {
    namespace.staticDeclarationsRecursive(coercedInputs, engineFunctions) match {
      case Success(d) => d.successNel
      case Failure(e) => s"Workflow has invalid declarations: ${e.getMessage}".failureNel
    }
  }

  private def validateNamespace(source: WdlSource): ErrorOr[NamespaceWithWorkflow] = {
    try {
      NamespaceWithWorkflow.load(source).successNel
    } catch {
      case e: Exception => s"Unable to load namespace from workflow: ${e.getMessage}".failureNel
    }
  }

  private def validateRawInputs(json: WdlJson): ErrorOr[Map[String, JsValue]] = {
    Try(json.parseJson) match {
      case Success(JsObject(inputs)) => inputs.successNel
      case Failure(reason: Throwable) => s"Workflow contains invalid inputs JSON: ${reason.getMessage}".failureNel
      case _ => s"Workflow inputs JSON cannot be parsed to JsObject: $json".failureNel
    }
  }

  private def validateCoercedInputs(rawInputs: Map[String, JsValue],
                                    namespace: NamespaceWithWorkflow): ErrorOr[WorkflowCoercedInputs] = {
    namespace.coerceRawInputs(rawInputs) match {
      case Success(r) => r.successNel
      case Failure(e: ThrowableWithErrors) => scalaz.Failure(e.errors)
      case Failure(e) => e.getMessage.failureNel
    }
  }

  private def validateWorkflowOptions(workflowOptions: WdlJson): ErrorOr[WorkflowOptions] = {
    def validateBackendOptions(backend: Backend, workflowOpt: WorkflowOptions): ErrorOr[WorkflowOptions] = {
      try {
        backend.assertWorkflowOptions(workflowOpt)
        workflowOpt.successNel
      } catch {
        case e: Exception => s"Workflow has invalid options for backend ${backend.backendType}: ${e.getMessage}".failureNel
      }
    }

    WorkflowOptions.fromJsonString(workflowOptions) match {
      case Success(options) =>
        val defaultBackend = CromwellBackend.getDefaultBackendForWorkflowAfterConsideringWorkflowOptions(options)
        // TODO: Also need to validate options against all assigned backends? Or all allowed backends?
        validateBackendOptions(defaultBackend, options)
      case Failure(e) => s"Workflow contains invalid options JSON: ${e.getMessage}".failureNel
    }
  }

  private def validateWorkflowFailureMode(workflowOptions: WorkflowOptions, conf: Config): ErrorOr[WorkflowFailureMode] = {
    val modeString: Try[String] = workflowOptions.get("workflowFailureMode") match {
      case Success(x) => Success(x)
      case Failure(_: OptionNotFoundException) => Success(conf.getStringOption("workflow-failure-mode") getOrElse DefaultWorkflowFailureMode)
      case Failure(t) => Failure(t)
    }

    modeString flatMap WorkflowFailureMode.tryParse match {
        case Success(mode) => mode.successNel
        case Failure(t) => t.getMessage.failureNel
    }
  }

  // TODO: With PBE, this should be defined in the backend.
  // TODO: Add CromwellRuntimeAttributes as a dependency for this actor (arg in ctor) in case is not moved to specific backend when PBE is merged.
  private def validateRuntimeAttributes(backendAssignments: Map[FullyQualifiedName, Backend], namespaceWithWorkflow: NamespaceWithWorkflow): ErrorOr[Seq[Set[String]]] = {
    TryUtil.sequence(namespaceWithWorkflow.workflow.calls map {
      call =>
        val backend = backendAssignments(call.fullyQualifiedName)
        CromwellRuntimeAttributes.validateKeys(call.task.runtimeAttributes.attrs.keySet, backend.backendType)
    }) match {
      case Success(validatedRuntimeAttrs) => validatedRuntimeAttrs.successNel
      case Failure(reason) => "Failed to validate runtime attributes.".failureNel
    }
  }
}
