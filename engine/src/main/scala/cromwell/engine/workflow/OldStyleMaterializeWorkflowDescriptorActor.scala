package cromwell.engine.workflow

import java.nio.file.Paths

import akka.actor.{Actor, Props}
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.wdl.OldWorkflowEngineFunctions
import cromwell.core.{ErrorOr, OptionNotFoundException, WorkflowId, WorkflowOptions, _}
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.util.TryUtil
import lenthall.config.ScalaConfig.EnhancedScalaConfig
import spray.json.{JsObject, _}
import wdl4s._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._
import scalaz.Validation
import scalaz.Validation.FlatMap._

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleMaterializeWorkflowDescriptorActor {
  def props(): Props = Props(new OldStyleMaterializeWorkflowDescriptorActor)

  sealed trait MaterializeWorkflowDescriptorActorMessage
  case class MaterializeWorkflow(id: WorkflowId,
                                 workflowSourceFiles: WorkflowSourceFiles,
                                 conf: Config = ConfigFactory.load) extends MaterializeWorkflowDescriptorActorMessage
  sealed trait MaterializationResult extends MaterializeWorkflowDescriptorActorMessage
  case class MaterializeWorkflowDescriptorSuccess(workflowDescriptor: OldStyleWorkflowDescriptor) extends MaterializationResult
  case class MaterializeWorkflowDescriptorFailure(reason: Throwable) extends Exception with MaterializationResult

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

  def workflowLogOptions(conf: Config): Option[OldStyleWorkflowLogOptions] = {
    for {
      workflowConfig <- conf.getConfigOption("workflow-options")
      dir <- workflowConfig.getStringOption("workflow-log-dir") if !dir.isEmpty
      temporary <- workflowConfig.getBooleanOption("workflow-log-temporary") orElse Option(true)
    } yield OldStyleWorkflowLogOptions(Paths.get(dir), temporary)
  }
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleMaterializeWorkflowDescriptorActor() extends Actor with LazyLogging {

  import OldStyleMaterializeWorkflowDescriptorActor._

  override def receive = {
    case MaterializeWorkflow(workflowId, workflowSourceFiles, conf) =>
      val backend = CromwellBackends.getBackendFromOptions(workflowSourceFiles.workflowOptionsJson)
      buildWorkflowDescriptor(workflowId, workflowSourceFiles, backend, conf) match {
        case scalaz.Success(descriptor) => sender() ! MaterializeWorkflowDescriptorSuccess(descriptor)
        case scalaz.Failure(error) => sender() ! MaterializeWorkflowDescriptorFailure(
          new IllegalArgumentException() with ThrowableWithErrors {
            val message = s"Workflow input processing failed."
            val errors = error
          })
      }
    case unknownMsg => logger.error(s"${this.getClass.getName} received an unknown message: $unknownMsg")
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFiles,
                                      backend: OldStyleBackend,
                                      conf: Config): ErrorOr[OldStyleWorkflowDescriptor] = {

    def buildWorkflowDescriptor(namespace: NamespaceWithWorkflow,
                                workflowOptions: WorkflowOptions,
                                rawInputs: Map[String, JsValue],
                                workflowFailureMode: WorkflowFailureMode): ErrorOr[OldStyleWorkflowDescriptor] = {
      validateCoercedInputs(rawInputs, namespace) flatMap { coercedInputs =>
        val workflowRootPath = backend.buildWorkflowRootPath(backend.rootPath(workflowOptions), namespace.workflow.unqualifiedName, id)
        val wfContext = new OldWorkflowContext(workflowRootPath)
        // PBE this is some disgustingness around the fact that creating the filesystems is the first
        // time we actually exercise the Google auth on the GCS filesystem, and if that auth is bad the
        // underlying code will throw.  This needs to be cleaned up as part of the explicit backend
        // initialization.  The auth used for 'genomics' also needs to be validated, which isn't addressed
        // by this hack.
        // More thinking is required about the relationship between filesystems and backends in general.
        val filesystemValidation = Validation.fromTryCatchNonFatal(backend.fileSystems(workflowOptions)).leftMap(_.getMessage).toValidationNel
        filesystemValidation flatMap { fileSystems =>
          val engineFunctions = backend.engineFunctions(fileSystems, wfContext)

          validateDeclarations(namespace, workflowOptions, coercedInputs, engineFunctions) flatMap { declarations =>
            OldStyleWorkflowDescriptor(id, sourceFiles, workflowOptions, workflowLogOptions(conf), rawInputs, namespace, coercedInputs, declarations, backend,
              configCallCaching(conf), lookupDockerHash(conf), workflowFailureMode, wfContext, fileSystems).successNel
          }
        }
      }
    }

    val namespaceValidation = validateNamespace(sourceFiles.wdlSource)
    val workflowOptionsValidation = validateWorkflowOptions(backend, sourceFiles.workflowOptionsJson)
    (namespaceValidation |@| workflowOptionsValidation) {
      (_, _)
    } flatMap { case (namespace, workflowOptions) =>
      val runtimeAttributes = validateRuntimeAttributes(namespace, backend.backendType)
      val rawInputsValidation = validateRawInputs(sourceFiles.inputsJson)
      val failureModeValidation = validateWorkflowFailureMode(workflowOptions, conf)
      (runtimeAttributes |@| rawInputsValidation |@| failureModeValidation) {
        (_, _, _)
      } flatMap { case (_, rawInputs, failureMode) =>
        buildWorkflowDescriptor(namespace, workflowOptions, rawInputs, failureMode)
      }
    }
  }

  private def validateDeclarations(namespace: NamespaceWithWorkflow,
                                   options: WorkflowOptions,
                                   coercedInputs: WorkflowCoercedInputs,
                                   engineFunctions: OldWorkflowEngineFunctions): ErrorOr[WorkflowCoercedInputs] = {
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

  // TODO: With PBE, this should be defined in the backend.
  private def validateWorkflowOptions(backend: OldStyleBackend, workflowOptions: WdlJson): ErrorOr[WorkflowOptions] = {
    def validateBackendOptions(backend: OldStyleBackend, workflowOpt: WorkflowOptions): ErrorOr[WorkflowOptions] = {
      try {
        backend.assertWorkflowOptions(workflowOpt)
        workflowOpt.successNel
      } catch {
        case e: Exception => s"Workflow has invalid options for backend ${backend.backendType}: ${e.getMessage}".failureNel
      }
    }

    WorkflowOptions.fromJsonString(workflowOptions) match {
      case Success(o) => validateBackendOptions(backend, o)
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
  // TODO: Add CromwellRuntimeAttributes as a dependency for this actor (arg in ctor) in case is not moved to specific
  // backend when PBE is merged.
  private def validateRuntimeAttributes(namespaceWithWorkflow: NamespaceWithWorkflow, backendType: BackendType): ErrorOr[Seq[Set[String]]] = {
    TryUtil.sequence(namespaceWithWorkflow.workflow.calls map {
      call => CromwellRuntimeAttributes.validateKeys(call.task.runtimeAttributes.attrs.keySet, backendType)
    }) match {
      case Success(validatedRuntimeAttrs) => validatedRuntimeAttrs.successNel
      case Failure(reason) => "Failed to validate runtime attributes.".failureNel
    }
  }
}
