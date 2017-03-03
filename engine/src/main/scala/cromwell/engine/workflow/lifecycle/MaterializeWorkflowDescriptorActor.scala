package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import cats.data.NonEmptyList
import cats.data.Validated._
import cats.instances.list._
import cats.instances.vector._
import cats.syntax.cartesian._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowOptions.{ReadFromCache, WorkflowOption, WriteToCache}
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.core.labels.{Label, Labels}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder}
import cromwell.engine._
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorActorData, MaterializeWorkflowDescriptorActorState}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import net.ceedubs.ficus.Ficus._
import spray.json._
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values.{WdlString, WdlValue}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object MaterializeWorkflowDescriptorActor {

  val RuntimeBackendKey: String = "backend"

  // This is a def so that the 'get' is only used when needed. And when it's needed, if the get fails
  // then initialization hasn't happened as we expected. As an indication that this is ok, previously
  // we might have called CromwellBackends.evaluateIfInitialized() which would have thrown a similar
  // exception if not initialized yet.
  def cromwellBackends = CromwellBackends.instance.get

  def props(serviceRegistryActor: ActorRef, workflowId: WorkflowId, cromwellBackends: => CromwellBackends = cromwellBackends, importLocalFilesystem: Boolean): Props = {
    Props(new MaterializeWorkflowDescriptorActor(serviceRegistryActor, workflowId, cromwellBackends, importLocalFilesystem)).withDispatcher(EngineDispatcher)
  }

  /*
  Commands
   */
  sealed trait MaterializeWorkflowDescriptorActorMessage
  case class MaterializeWorkflowDescriptorCommand(workflowSourceFiles: WorkflowSourceFilesCollection,
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
  case class MaterializeWorkflowDescriptorActorData()

  private val DefaultWorkflowFailureMode = NoNewCalls.toString

  private[lifecycle] def validateCallCachingMode(workflowOptions: WorkflowOptions, conf: Config): ErrorOr[CallCachingMode] = {

    def readOptionalOption(option: WorkflowOption): ErrorOr[Boolean] = {
      workflowOptions.getBoolean(option.name) match {
        case Success(x) => x.validNel
        case Failure(_: OptionNotFoundException) => true.validNel
        case Failure(t) => t.getMessage.invalidNel
      }
    }

    val enabled = conf.as[Option[Boolean]]("call-caching.enabled").getOrElse(false)
    val invalidateBadCacheResults = conf.as[Option[Boolean]]("call-caching.invalidate-bad-cache-results").getOrElse(true)
    val callCachingOptions = CallCachingOptions(invalidateBadCacheResults)
    if (enabled) {
      val readFromCache = readOptionalOption(ReadFromCache)
      val writeToCache = readOptionalOption(WriteToCache)

      (readFromCache |@| writeToCache) map {
        case (false, false) => CallCachingOff
        case (true, false) => CallCachingActivity(ReadCache, callCachingOptions)
        case (false, true) => CallCachingActivity(WriteCache, callCachingOptions)
        case (true, true) => CallCachingActivity(ReadAndWriteCache, callCachingOptions)
      }
    }
    else {
      CallCachingOff.validNel
    }
  }
}

class MaterializeWorkflowDescriptorActor(serviceRegistryActor: ActorRef,
                                         val workflowIdForLogging: WorkflowId,
                                         cromwellBackends: => CromwellBackends,
                                         importLocalFilesystem: Boolean) extends LoggingFSM[MaterializeWorkflowDescriptorActorState, MaterializeWorkflowDescriptorActorData] with LazyLogging with WorkflowLogging {

  import MaterializeWorkflowDescriptorActor._

  val tag = self.path.name

  val iOExecutionContext = context.system.dispatchers.lookup("akka.dispatchers.io-dispatcher")

  startWith(ReadyToMaterializeState, MaterializeWorkflowDescriptorActorData())

  when(ReadyToMaterializeState) {
    case Event(MaterializeWorkflowDescriptorCommand(workflowSourceFiles, conf), _) =>
      buildWorkflowDescriptor(workflowIdForLogging, workflowSourceFiles, conf) match {
        case Valid(descriptor) =>
          sender() ! MaterializeWorkflowDescriptorSuccessResponse(descriptor)
          goto(MaterializationSuccessfulState)
        case Invalid(error) =>
          sender() ! MaterializeWorkflowDescriptorFailureResponse(
            new IllegalArgumentException with MessageAggregation {
              val exceptionContext = s"Workflow input processing failed"
              val errorMessages = error.toList
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
      workflowLogger.debug(s"transition from $oldState to $terminalState. Stopping self")
      context.stop(self)
    case fromState -> toState =>
      workflowLogger.debug(s"transitioning from $fromState to $toState")
  }

  whenUnhandled {
    case Event(EngineLifecycleActorAbortCommand, _) =>
      goto(MaterializationAbortedState)
    case unhandledMessage =>
      workflowLogger.warn(s"received an unhandled message $unhandledMessage in state $stateName")
      stay
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFilesCollection,
                                      conf: Config): ErrorOr[EngineWorkflowDescriptor] = {
    val namespaceValidation = validateNamespace(sourceFiles)
    val workflowOptionsValidation = validateWorkflowOptions(sourceFiles.workflowOptionsJson)
    val labelsValidation = validateLabels(sourceFiles.labelsJson)
    (namespaceValidation |@| workflowOptionsValidation |@| labelsValidation) map {
      (_, _, _)
    } flatMap { case (namespace, workflowOptions, labels) =>
      pushWfNameMetadataService(namespace.workflow.unqualifiedName)
      val pathBuilders = EngineFilesystems(context.system).pathBuildersForWorkflow(workflowOptions)
      buildWorkflowDescriptor(id, sourceFiles, namespace, workflowOptions, labels, conf, pathBuilders)
    }
  }


  private def pushWfNameMetadataService(name: String): Unit = {
    // Workflow name:
    val nameEvent = MetadataEvent(MetadataKey(workflowIdForLogging, None, WorkflowMetadataKeys.Name), MetadataValue(name))

    serviceRegistryActor ! PutMetadataAction(nameEvent)
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFilesCollection,
                                      namespace: WdlNamespaceWithWorkflow,
                                      workflowOptions: WorkflowOptions,
                                      labels: Labels,
                                      conf: Config,
                                      pathBuilders: List[PathBuilder]): ErrorOr[EngineWorkflowDescriptor] = {
    val defaultBackendName = conf.as[Option[String]]("backend.default")
    val rawInputsValidation = validateRawInputs(sourceFiles.inputsJson)

    val failureModeValidation = validateWorkflowFailureMode(workflowOptions, conf)
    val backendAssignmentsValidation = validateBackendAssignments(namespace.taskCalls, workflowOptions, defaultBackendName)
    val callCachingModeValidation = validateCallCachingMode(workflowOptions, conf)

    (rawInputsValidation |@| failureModeValidation |@| backendAssignmentsValidation |@| callCachingModeValidation) map {
      (_, _, _, _)
    } flatMap { case (rawInputs, failureMode, backendAssignments, callCachingMode) =>
      buildWorkflowDescriptor(id, namespace, rawInputs, backendAssignments, workflowOptions, labels, failureMode, pathBuilders, callCachingMode)
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      namespace: WdlNamespaceWithWorkflow,
                                      rawInputs: Map[String, JsValue],
                                      backendAssignments: Map[TaskCall, String],
                                      workflowOptions: WorkflowOptions,
                                      labels: Labels,
                                      failureMode: WorkflowFailureMode,
                                      pathBuilders: List[PathBuilder],
                                      callCachingMode: CallCachingMode): ErrorOr[EngineWorkflowDescriptor] = {

    def checkTypes(inputs: Map[FullyQualifiedName, WdlValue]): ErrorOr[Map[FullyQualifiedName, WdlValue]] = {
      val allDeclarations = namespace.workflow.declarations ++ namespace.workflow.calls.flatMap(_.declarations)
      val list: List[ErrorOr[(FullyQualifiedName, WdlValue)]] = inputs.map({ case (k, v) =>
        allDeclarations.find(_.fullyQualifiedName == k) match {
          case Some(decl) if decl.wdlType.coerceRawValue(v).isFailure =>
            s"Invalid right-side type of '$k'.  Expecting ${decl.wdlType.toWdlString}, got ${v.wdlType.toWdlString}".invalidNel
          case _ => (k, v).validNel[String]
        }
      }).toList

      val validatedInputs: ErrorOr[List[(FullyQualifiedName, WdlValue)]] = list.sequence[ErrorOr, (FullyQualifiedName, WdlValue)]
      validatedInputs.map(_.toMap)
    }

    for {
      coercedInputs <- validateCoercedInputs(rawInputs, namespace)
      _ = pushWfInputsToMetadataService(coercedInputs)
      evaluatedWorkflowsDeclarations <- validateDeclarations(namespace, workflowOptions, coercedInputs, pathBuilders)
      declarationsAndInputs <- checkTypes(evaluatedWorkflowsDeclarations ++ coercedInputs)
      backendDescriptor = BackendWorkflowDescriptor(id, namespace.workflow, declarationsAndInputs, workflowOptions, labels)
    } yield EngineWorkflowDescriptor(namespace, backendDescriptor, backendAssignments, failureMode, pathBuilders, callCachingMode)
  }

  private def pushWfInputsToMetadataService(workflowInputs: WorkflowCoercedInputs): Unit = {
    // Inputs
    val inputEvents = workflowInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(MetadataKey(workflowIdForLogging, None,WorkflowMetadataKeys.Inputs)))
      case inputs =>
        inputs flatMap { case (inputName, wdlValue) =>
          wdlValueToMetadataEvents(MetadataKey(workflowIdForLogging, None, s"${WorkflowMetadataKeys.Inputs}:$inputName"), wdlValue)
        }
    }

    serviceRegistryActor ! PutMetadataAction(inputEvents)
  }

  private def validateBackendAssignments(calls: Set[TaskCall], workflowOptions: WorkflowOptions, defaultBackendName: Option[String]): ErrorOr[Map[TaskCall, String]] = {
    val callToBackendMap = Try {
      calls map { call =>
        val backendPriorities = Seq(
          workflowOptions.get(RuntimeBackendKey).toOption,
          assignBackendUsingRuntimeAttrs(call),
          defaultBackendName
        )

        backendPriorities.flatten.headOption match {
          case Some(backendName) if cromwellBackends.isValidBackendName(backendName) => call -> backendName
          case Some(backendName) => throw new Exception(s"Backend for call ${call.fullyQualifiedName} ('$backendName') not registered in configuration file")
          case None => throw new Exception(s"No backend could be found for call ${call.fullyQualifiedName}")
        }
      } toMap
    }

    callToBackendMap match {
      case Success(backendMap) =>
        val backendMapAsString = backendMap.map({case (k, v) => s"${k.fullyQualifiedName} -> $v"}).mkString(", ")
        workflowLogger.info(s"Call-to-Backend assignments: $backendMapAsString")
        backendMap.validNel
      case Failure(t) => t.getMessage.invalidNel
    }
  }

  /**
    * Map a call to a backend name depending on the runtime attribute key
    */
  private def assignBackendUsingRuntimeAttrs(call: TaskCall): Option[String] = {
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

  private def validateDeclarations(namespace: WdlNamespaceWithWorkflow,
                                   options: WorkflowOptions,
                                   coercedInputs: WorkflowCoercedInputs,
                                   pathBuilders: List[PathBuilder]): ErrorOr[WorkflowCoercedInputs] = {
    namespace.staticDeclarationsRecursive(coercedInputs, new WdlFunctions(pathBuilders)) match {
      case Success(d) => d.validNel
      case Failure(e) => s"Workflow has invalid declarations: ${e.getMessage}".invalidNel
    }
  }

  private def validateImportsDirectory(zipContents: Array[Byte]): ErrorOr[Path] = {

    def makeZipFile(contents: Array[Byte]): Try[Path] = Try {
      DefaultPathBuilder.createTempFile("", ".zip").write(contents)(OpenOptions.default)
    }

    def unZipFile(f: Path) = Try {
      val unzippedFile = f.unzip()
      val unzippedFileContents = unzippedFile.list.toSeq.head
      if (unzippedFileContents.isDirectory) unzippedFileContents else unzippedFile
    }

    val importsFile = for {
      zipFile <- makeZipFile(zipContents)
      unzipped <- unZipFile(zipFile)
      _ <- Try(zipFile.delete(swallowIOExceptions = true))
    } yield unzipped

    importsFile match {
      case Success(unzippedDirectory: Path) => unzippedDirectory.validNel
      case Failure(t) => t.getMessage.invalidNel
    }
  }

  private def validateNamespaceWithImports(w: WorkflowSourceFilesWithDependenciesZip): ErrorOr[WdlNamespaceWithWorkflow] = {
    def getMetadatae(importsDir: Path, prefix: String = ""): List[(String, Path)] = {
      importsDir.children.toList flatMap {
        case f: Path if f.isDirectory => getMetadatae(f, prefix + f.name + "/")
        case f: Path if f.name.endsWith(".wdl") => List((prefix + f.name, f))
        case _ => List.empty
      }
    }

    def writeMetadatae(importsDir: Path) = {
      val wfImportEvents = getMetadatae(importsDir) map { case (name: String, f: Path) =>
        val contents = f.lines.mkString(System.lineSeparator())
        MetadataEvent(MetadataKey(workflowIdForLogging, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Imports, name), MetadataValue(contents))
      }
      serviceRegistryActor ! PutMetadataAction(wfImportEvents)
    }

    def importsAsNamespace(importsDir: Path): ErrorOr[WdlNamespaceWithWorkflow] = {
      writeMetadatae(importsDir)
      val importsDirFile = better.files.File(importsDir.pathAsString) // For wdl4s better file compatibility
      val importResolvers: Seq[ImportResolver] = if (importLocalFilesystem) {
        List(WdlNamespace.directoryResolver(importsDirFile), WdlNamespace.fileResolver)
      } else {
        List(WdlNamespace.directoryResolver(importsDirFile))
      }
      val results = WdlNamespaceWithWorkflow.load(w.wdlSource, importResolvers)
      importsDir.delete(swallowIOExceptions = true)
      results match {
        case Success(ns) => ns.validNel
        case Failure(f) => f.getMessage.invalidNel
      }
    }

    validateImportsDirectory(w.importsZip) flatMap importsAsNamespace
  }

  private def validateNamespace(source: WorkflowSourceFilesCollection): ErrorOr[WdlNamespaceWithWorkflow] = {
    try {
      source match {
        case w: WorkflowSourceFilesWithDependenciesZip => validateNamespaceWithImports(w)
        case w: WorkflowSourceFilesWithoutImports =>
          val importResolvers: Seq[ImportResolver] = if (importLocalFilesystem) {
            List(WdlNamespace.fileResolver)
          } else {
            List.empty
          }
          // This .get is ok because we're already in a try/catch.
          WdlNamespaceWithWorkflow.load(w.wdlSource, importResolvers).get.validNel
      }
    } catch {
      case e: Exception => s"Unable to load namespace from workflow: ${e.getMessage}".invalidNel
    }
  }

  private def validateRawInputs(json: WdlJson): ErrorOr[Map[String, JsValue]] = {
    Try(json.parseJson) match {
      case Success(JsObject(inputs)) => inputs.validNel
      case Failure(reason: Throwable) => s"Workflow contains invalid inputs JSON: ${reason.getMessage}".invalidNel
      case _ => s"Workflow inputs JSON cannot be parsed to JsObject: $json".invalidNel
    }
  }

  private def validateLabels(json: WdlJson): ErrorOr[Labels] = {

    def toLabels(inputs: Map[String, JsValue]): ErrorOr[Labels] = {
      val vectorOfValidatedLabel: Vector[ErrorOr[Label]] = inputs.toVector map {
        case (key, JsString(s)) => Label.validateLabel(key, s)
        case (key, other) => s"Invalid label '$key: $other': Labels must be strings, and must match the regex ${Label.LabelRegexPattern}".invalidNel
      }

      vectorOfValidatedLabel.sequence[ErrorOr, Label] map { validatedVectorofLabel => Labels(validatedVectorofLabel) }
    }

    Try(json.parseJson) match {
      case Success(JsObject(inputs)) => toLabels(inputs)
      case Failure(reason: Throwable) => s"Workflow contains invalid labels JSON: ${reason.getMessage}".invalidNel
      case _ => """Invalid workflow labels JSON. Expected a JsObject of "labelKey": "labelValue" values.""".invalidNel
    }
  }

  private def validateCoercedInputs(rawInputs: Map[String, JsValue],
                                    namespace: WdlNamespaceWithWorkflow): ErrorOr[WorkflowCoercedInputs] = {
    namespace.coerceRawInputs(rawInputs) match {
      case Success(r) => r.validNel
      case Failure(e: MessageAggregation) if e.errorMessages.nonEmpty => Invalid(NonEmptyList.fromListUnsafe(e.errorMessages.toList))
      case Failure(e) => e.getMessage.invalidNel
    }
  }

  private def validateWorkflowOptions(workflowOptions: WdlJson): ErrorOr[WorkflowOptions] = {
    WorkflowOptions.fromJsonString(workflowOptions) match {
      case Success(opts) => opts.validNel
      case Failure(e) => s"Workflow contains invalid options JSON: ${e.getMessage}".invalidNel
    }
  }

  private def validateWorkflowFailureMode(workflowOptions: WorkflowOptions, conf: Config): ErrorOr[WorkflowFailureMode] = {
    val modeString: Try[String] = workflowOptions.get(WorkflowOptions.WorkflowFailureMode) match {
      case Success(x) => Success(x)
      case Failure(_: OptionNotFoundException) => Success(conf.as[Option[String]]("workflow-options.workflow-failure-mode") getOrElse DefaultWorkflowFailureMode)
      case Failure(t) => Failure(t)
    }

    modeString flatMap WorkflowFailureMode.tryParse match {
        case Success(mode) => mode.validNel
        case Failure(t) => t.getMessage.invalidNel
    }
  }
}
