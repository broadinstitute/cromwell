package cromwell.engine.workflow.lifecycle

import java.nio.file.Paths

import akka.actor.{FSM, LoggingFSM, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core._
import cromwell.engine._
import cromwell.engine.backend.{CromwellBackends, OldStyleWorkflowLogOptions}
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorActorState, ShadowMaterializeWorkflowDescriptorActorData}
import lenthall.config.ScalaConfig.EnhancedScalaConfig
import spray.json.{JsObject, _}
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values.{WdlString, WdlValue}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._

object MaterializeWorkflowDescriptorActor {

  val RuntimeBackendKey: String = "backend"

  def props(): Props = Props(new MaterializeWorkflowDescriptorActor)

  /*
  Commands
   */
  sealed trait MaterializeWorkflowDescriptorActorMessage
  case class MaterializeWorkflowDescriptorCommand(id: WorkflowId,
                                                  workflowSourceFiles: WorkflowSourceFiles,
                                                  conf: Config) extends MaterializeWorkflowDescriptorActorMessage
  case object MaterializeWorkflowDescriptorAbortCommand

  /*
  Responses
   */
  sealed trait WorkflowDescriptorMaterializationResult extends MaterializeWorkflowDescriptorActorMessage
  case class MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor: EngineWorkflowDescriptor) extends WorkflowDescriptorMaterializationResult
  case class MaterializeWorkflowDescriptorFailureResponse(reason: Throwable) extends Exception with WorkflowDescriptorMaterializationResult

  /*
  States
   */
  sealed trait MaterializeWorkflowDescriptorActorState { def terminal = false }
  sealed trait MaterializeWorkflowDescriptorActorTerminalState extends MaterializeWorkflowDescriptorActorState {
    override val terminal = true
  }
  case object ReadyToMaterializeState extends MaterializeWorkflowDescriptorActorState
  case object MaterializationSuccessfulState extends MaterializeWorkflowDescriptorActorTerminalState
  case object MaterializationFailedState extends MaterializeWorkflowDescriptorActorTerminalState
  case object MaterializationAbortedState extends MaterializeWorkflowDescriptorActorTerminalState

  /*
  Data
   */
  case class ShadowMaterializeWorkflowDescriptorActorData()

  import lenthall.config.ScalaConfig._

  private val DefaultWorkflowFailureMode = NoNewCalls.toString

  def workflowLogOptions(conf: Config): Option[OldStyleWorkflowLogOptions] = {
    for {
      workflowConfig <- conf.getConfigOption("workflow-options")
      dir <- workflowConfig.getStringOption("workflow-log-dir") if !dir.isEmpty
      temporary <- workflowConfig.getBooleanOption("workflow-log-temporary") orElse Option(true)
    } yield OldStyleWorkflowLogOptions(Paths.get(dir), temporary)
  }
}

class MaterializeWorkflowDescriptorActor() extends LoggingFSM[MaterializeWorkflowDescriptorActorState, ShadowMaterializeWorkflowDescriptorActorData] with LazyLogging {

  import MaterializeWorkflowDescriptorActor._

  val tag = self.path.name

  startWith(ReadyToMaterializeState, ShadowMaterializeWorkflowDescriptorActorData())

  when(ReadyToMaterializeState) {
    case Event(MaterializeWorkflowDescriptorCommand(workflowId, workflowSourceFiles, conf), _) =>
      buildWorkflowDescriptor(workflowId, workflowSourceFiles, conf) match {
        case scalaz.Success(descriptor) =>
          sender() ! MaterializeWorkflowDescriptorSuccessResponse(descriptor)
          goto(MaterializationSuccessfulState)
        case scalaz.Failure(error) =>
          sender() ! MaterializeWorkflowDescriptorFailureResponse(new IllegalArgumentException() with ThrowableWithErrors {
            val message = s"Workflow input processing failed."
            val errors = error
          })
          goto(MaterializationFailedState)
      }
    case Event(MaterializeWorkflowDescriptorAbortCommand, _) =>
      goto(MaterializationAbortedState)
  }

  // Let these fall through to the whenUnhandled handler:
  when(MaterializationSuccessfulState) { FSM.NullFunction }
  when(MaterializationFailedState) { FSM.NullFunction }
  when(MaterializationAbortedState) { FSM.NullFunction }

  onTransition {
    case oldState -> terminalState if terminalState.terminal =>
      log.info(s"$tag transition from $oldState to $terminalState: shutting down")
      context.stop(self)
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState")
  }

  whenUnhandled {
    case Event(EngineLifecycleActorAbortCommand, _) =>
      goto(MaterializationAbortedState)
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message $unhandledMessage in state $stateName")
      stay
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFiles,
                                      conf: Config): ErrorOr[EngineWorkflowDescriptor] = {
    val namespaceValidation = validateNamespace(sourceFiles.wdlSource)
    val workflowOptionsValidation = validateWorkflowOptions(sourceFiles.workflowOptionsJson)
    (namespaceValidation |@| workflowOptionsValidation) {
      (_, _)
    } flatMap { case (namespace, workflowOptions) =>
      buildWorkflowDescriptor(id, sourceFiles, namespace, workflowOptions, conf)
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFiles,
                                      namespace: NamespaceWithWorkflow,
                                      workflowOptions: WorkflowOptions,
                                      conf: Config): ErrorOr[EngineWorkflowDescriptor] = {
    val rawInputsValidation = validateRawInputs(sourceFiles.inputsJson)
    val failureModeValidation = validateWorkflowFailureMode(workflowOptions, conf)
    val backendAssignmentsValidation = validateBackendAssignments(namespace.workflow.calls, workflowOptions)
    (rawInputsValidation |@| failureModeValidation |@| backendAssignmentsValidation ) {
      (_, _, _)
    } flatMap { case (rawInputs, failureMode, backendAssignments) =>
      buildWorkflowDescriptor(id, namespace, rawInputs, backendAssignments, workflowOptions, failureMode)
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      namespace: NamespaceWithWorkflow,
                                      rawInputs: Map[String, JsValue],
                                      backendAssignments: Map[Call, String],
                                      workflowOptions: WorkflowOptions,
                                      failureMode: WorkflowFailureMode): ErrorOr[EngineWorkflowDescriptor] = {
    for {
      coercedInputs <- validateCoercedInputs(rawInputs, namespace)
      declarations <- validateDeclarations(namespace, workflowOptions, coercedInputs)
      declarationsAndInputs = declarations ++ coercedInputs
      backendDescriptor = BackendWorkflowDescriptor(id, namespace, declarationsAndInputs, workflowOptions)
    } yield EngineWorkflowDescriptor(backendDescriptor, declarations, backendAssignments, failureMode)
  }

  private def validateBackendAssignments(calls: Seq[Call], workflowOptions: WorkflowOptions): ErrorOr[Map[Call, String]] = {
    val callToBackendMap = Try {
      calls map { call =>
        val backendPriorities = Seq(
          workflowOptions.get(RuntimeBackendKey).toOption,
          assignBackendUsingRuntimeAttrs(call),
          Option(CromwellBackends.shadowDefaultBackend)
        )

        backendPriorities.flatten.headOption match {
          case Some(backendName) if CromwellBackends.isValidBackendName(backendName) => call -> backendName
          case Some(backendName) => throw new Exception(s"Invalid backend for call ${call.fullyQualifiedName}: $backendName")
          case None => throw new Exception(s"No backend could be found for call ${call.fullyQualifiedName}")
        }
      } toMap
    }

    callToBackendMap match {
      case Success(backendMap) =>
        val backendMapAsString = backendMap.map({case (k, v) => s"${k.fullyQualifiedName} -> $v"}).mkString(", ")
        log.info(s"$tag: Call-to-Backend assignments: $backendMapAsString")
        backendMap.successNel
      case Failure(t) => t.getMessage.failureNel
    }
  }

  /**
    * Map a call to a backend name depending on the runtime attribute key
    */
  private def assignBackendUsingRuntimeAttrs(call: Call): Option[String] = {
    val runtimeAttributesMap = call.task.runtimeAttributes.attrs
    runtimeAttributesMap.get(RuntimeBackendKey) map { wdlExpr => evaluateBackendNameExpression(call.fullyQualifiedName, wdlExpr) }
  }

  private def evaluateBackendNameExpression(callName: String, backendNameAsExp: WdlExpression): String = {
    backendNameAsExp.evaluate(NoLookup, NoFunctions) match {
      case Success(runtimeString: WdlString) => runtimeString.valueString
      case Success(x: WdlValue) =>
        throw new Exception(s"Non-string values are not currently supported for backends! Cannot use backend '${x.valueString}' to backend to Call: $callName")
      case Failure(error) =>
        throw new Exception(s"Dynamic backends are not currently supported! Cannot assign backend '${backendNameAsExp.valueString}' for Call: $callName", error)
    }
  }

  private def validateDeclarations(namespace: NamespaceWithWorkflow,
                                   options: WorkflowOptions,
                                   coercedInputs: WorkflowCoercedInputs): ErrorOr[WorkflowCoercedInputs] = {
    // TODO: Need to create engine-only engine functions!
    namespace.staticWorkflowDeclarationsRecursive(coercedInputs, new WdlFunctions(options)) match {
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
    WorkflowOptions.fromJsonString(workflowOptions) match {
      case Success(opts) => opts.successNel
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
}
