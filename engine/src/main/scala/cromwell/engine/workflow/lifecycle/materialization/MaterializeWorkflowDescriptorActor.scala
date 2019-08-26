package cromwell.engine.workflow.lifecycle.materialization

import akka.actor.{ActorRef, FSM, LoggingFSM, Props, Status}
import akka.pattern.pipe
import cats.data.EitherT._
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import common.exception.{AggregatedMessageException, MessageAggregation}
import common.validation.Checked._
import common.validation.ErrorOr._
import common.validation.IOChecked._
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowOptions.{ReadFromCache, WorkflowOption, WriteToCache}
import cromwell.core.callcaching._
import cromwell.core.io.AsyncIo
import cromwell.core.labels.{Label, Labels}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.path.{PathBuilder, PathBuilderFactory}
import cromwell.core._
import cromwell.engine._
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.WorkflowProcessingEventPublishing._
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.languages.util.ImportResolver._
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import net.ceedubs.ficus.Ficus._
import spray.json._
import wom.core.WorkflowSource
import wom.expression.{NoIoFunctionSet, WomExpression}
import wom.graph.CommandCallNode
import wom.graph.GraphNodePort.OutputPort
import wom.runtime.WomOutputRuntimeExtractor
import wom.values.{WomString, WomValue}

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object MaterializeWorkflowDescriptorActor {

  val RuntimeBackendKey: String = "backend"

  // This is a def so that the 'get' is only used when needed. And when it's needed, if the get fails
  // then initialization hasn't happened as we expected. As an indication that this is ok, previously
  // we might have called CromwellBackends.evaluateIfInitialized() which would have thrown a similar
  // exception if not initialized yet.
  private def cromwellBackends = CromwellBackends.instance.get

  def props(serviceRegistryActor: ActorRef, workflowId: WorkflowId, cromwellBackends: => CromwellBackends = cromwellBackends,
            importLocalFilesystem: Boolean, ioActorProxy: ActorRef, hogGroup: HogGroup): Props = {
    Props(new MaterializeWorkflowDescriptorActor(serviceRegistryActor, workflowId, cromwellBackends, importLocalFilesystem, ioActorProxy, hogGroup)).withDispatcher(EngineDispatcher)
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
  case object MaterializingState extends MaterializeWorkflowDescriptorActorState
  case object MaterializationSuccessfulState extends MaterializeWorkflowDescriptorActorTerminalState
  case object MaterializationFailedState extends MaterializeWorkflowDescriptorActorTerminalState
  case object MaterializationAbortedState extends MaterializeWorkflowDescriptorActorTerminalState

  private val DefaultWorkflowFailureMode = NoNewCalls.toString

  private[lifecycle] def validateCallCachingMode(workflowOptions: WorkflowOptions, conf: Config): ErrorOr[CallCachingMode] = {

    def readOptionalOption(option: WorkflowOption): ErrorOr[Boolean] = {
      workflowOptions.getBoolean(option.name) match {
        case Success(x) => x.validNel
        case Failure(_: OptionNotFoundException) => true.validNel
        case Failure(t) => t.getMessage.invalidNel
      }
    }

    val callCachingConfig = conf.getConfig("call-caching")

    def errorOrCallCachingBoolean(path: String): ErrorOr[Boolean] = {
      import common.validation.Validation._
      validate(callCachingConfig.getBoolean(path))
    }

    val errorOrEnabled = errorOrCallCachingBoolean("enabled")
    if (errorOrEnabled.exists(_ == true)) {
      val readFromCache = readOptionalOption(ReadFromCache)
      val writeToCache = readOptionalOption(WriteToCache)

      def errorOrCallCachingMode(callCachingOptions: CallCachingOptions): ErrorOr[CallCachingMode] = {
        (readFromCache, writeToCache) mapN {
          case (false, false) => CallCachingOff
          case (true, false) => CallCachingActivity(ReadCache, callCachingOptions)
          case (false, true) => CallCachingActivity(WriteCache, callCachingOptions)
          case (true, true) => CallCachingActivity(ReadAndWriteCache, callCachingOptions)
        }
      }

      val errorOrMaybePrefixes = workflowOptions.getVectorOfStrings("call_cache_hit_path_prefixes")
      val errorOrInvalidateBadCacheResults = errorOrCallCachingBoolean("invalidate-bad-cache-results")
      val errorOrCallCachingOptions = (
        errorOrMaybePrefixes,
        errorOrInvalidateBadCacheResults,
      ) mapN {
        (
          maybePrefixes,
          invalidateBadCacheResults,
        ) =>
          CallCachingOptions(
            invalidateBadCacheResults,
            maybePrefixes,
          )
      }
      for {
        callCachingOptions <- errorOrCallCachingOptions
        mode <- errorOrCallCachingMode(callCachingOptions)
      } yield mode
    }
    else {
      errorOrEnabled.map(_ => CallCachingOff)
    }
  }
}

// TODO WOM: need to decide where to draw the line between language specific initialization and WOM
class MaterializeWorkflowDescriptorActor(serviceRegistryActor: ActorRef,
                                         workflowId: WorkflowId,
                                         cromwellBackends: => CromwellBackends,
                                         importLocalFilesystem: Boolean,
                                         ioActorProxy: ActorRef,
                                         hogGroup: HogGroup) extends LoggingFSM[MaterializeWorkflowDescriptorActorState, Unit] with StrictLogging with WorkflowLogging {

  import MaterializeWorkflowDescriptorActor._
  val tag = self.path.name

  override lazy val workflowIdForLogging = workflowId.toPossiblyNotRoot
  override lazy val rootWorkflowIdForLogging = workflowId.toRoot

  val iOExecutionContext = context.system.dispatchers.lookup("akka.dispatchers.io-dispatcher")
  implicit val ec = context.dispatcher

  protected val pathBuilderFactories: List[PathBuilderFactory] = EngineFilesystems.configuredPathBuilderFactories

  startWith(ReadyToMaterializeState, ())

  when(ReadyToMaterializeState) {
    case Event(MaterializeWorkflowDescriptorCommand(workflowSourceFiles, conf), _) =>
      val replyTo = sender()

      workflowOptionsAndPathBuilders(workflowSourceFiles) match {
        case (workflowOptions, pathBuilders) =>
          val futureDescriptor: Future[ErrorOr[EngineWorkflowDescriptor]] = pathBuilders flatMap { pb =>
            val engineIoFunctions = new EngineIoFunctions(pb, new AsyncIo(ioActorProxy, GcsBatchCommandBuilder), iOExecutionContext)
            buildWorkflowDescriptor(workflowId, workflowSourceFiles, conf, workflowOptions, pb, engineIoFunctions).
              value.
              unsafeToFuture().
              map(_.toValidated)
          }

          // Pipe the response to self, but make it look like it comes from the sender of the command
          // This way we can access it through sender() in the next state and don't have to store the value
          // of replyTo in the data
          pipe(futureDescriptor).to(self, replyTo)
          goto(MaterializingState)
      }
    case Event(MaterializeWorkflowDescriptorAbortCommand, _) =>
      goto(MaterializationAbortedState)
  }

  when(MaterializingState) {
    case Event(Valid(descriptor: EngineWorkflowDescriptor), _) =>
      sender() ! MaterializeWorkflowDescriptorSuccessResponse(descriptor)
      goto(MaterializationSuccessfulState)
    case Event(Invalid(error: NonEmptyList[String]@unchecked), _) =>
      workflowInitializationFailed(error, sender())
      goto(MaterializationFailedState)
    case Event(Status.Failure(failure), _) =>
      workflowInitializationFailed(NonEmptyList.of(failure.getMessage, failure.getStackTrace.map(_.toString):_*), sender())
      goto(MaterializationFailedState)
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

  private def workflowInitializationFailed(errors: NonEmptyList[String], replyTo: ActorRef) = {
    sender() ! MaterializeWorkflowDescriptorFailureResponse(
      new IllegalArgumentException with MessageAggregation {
        val exceptionContext = "Workflow input processing failed"
        val errorMessages = errors.toList
      })
  }

  private def workflowOptionsAndPathBuilders(sourceFiles: WorkflowSourceFilesCollection): (WorkflowOptions, Future[List[PathBuilder]]) = {
    sourceFiles.workflowOptions
    val pathBuilders = EngineFilesystems.pathBuildersForWorkflow(sourceFiles.workflowOptions, pathBuilderFactories)(context.system)
    (sourceFiles.workflowOptions, pathBuilders)
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFilesCollection,
                                      conf: Config,
                                      workflowOptions: WorkflowOptions,
                                      pathBuilders: List[PathBuilder],
                                      engineIoFunctions: EngineIoFunctions): IOChecked[EngineWorkflowDescriptor] = {

    def findFactory(workflowSource: WorkflowSource): ErrorOr[LanguageFactory] = {

      val factory = LanguageFactoryUtil.chooseFactory(workflowSource, sourceFiles)

      factory foreach { validFactory =>
        workflowLogger.info(s"Parsing workflow as ${validFactory.languageName} ${validFactory.languageVersionName}")
        pushLanguageToMetadata(validFactory.languageName, validFactory.languageVersionName)
      }
      
      factory
    }

    def buildValidatedNamespace(factory: LanguageFactory, workflowSource: WorkflowSource, importResolvers: List[ImportResolver]): IOChecked[ValidatedWomNamespace] = {
      factory.validateNamespace(
        sourceFiles,
        workflowSource,
        workflowOptions,
        importLocalFilesystem,
        workflowId,
        engineIoFunctions,
        importResolvers
      )
    }

    val localFilesystemResolvers =
      if (importLocalFilesystem) DirectoryResolver.localFilesystemResolvers(None)
      else List.empty

    val zippedResolverCheck: IOChecked[Option[DirectoryResolver]] = fromEither[IO](sourceFiles.importsZipFileOption match {
      case None => None.validNelCheck
      case Some(zipContent) => zippedImportResolver(zipContent, workflowId).toEither.map(Option.apply)
    })

    val labels = convertJsonToLabels(sourceFiles.labelsJson)

    for {
      _ <- publishLabelsToMetadata(id, labels.asMap, serviceRegistryActor)
      zippedImportResolver <- zippedResolverCheck
      importResolvers = zippedImportResolver.toList ++ localFilesystemResolvers :+ HttpResolver(None, Map.empty)
      sourceAndResolvers <- fromEither[IO](LanguageFactoryUtil.findWorkflowSource(sourceFiles.workflowSource, sourceFiles.workflowUrl, importResolvers))
      _ = if(sourceFiles.workflowUrl.isDefined) publishWorkflowSourceToMetadata(id, sourceAndResolvers._1)
      factory <- findFactory(sourceAndResolvers._1).toIOChecked
      outputRuntimeExtractor <- factory.womOutputRuntimeExtractor.toValidated.toIOChecked
      validatedNamespace <- buildValidatedNamespace(factory, sourceAndResolvers._1, sourceAndResolvers._2)
      closeResult = sourceAndResolvers._2.traverse(_.cleanupIfNecessary())
      _ = pushNamespaceMetadata(validatedNamespace)
      ewd <- buildWorkflowDescriptor(id, validatedNamespace, workflowOptions, labels, conf, pathBuilders, outputRuntimeExtractor).toIOChecked
    } yield {
      closeResult match {
        case Valid(_) => ()
        case Invalid(errorsList) =>
          logger.warn(s"Import resolver(s) closed with errors: [${errorsList.toList.mkString(", ")}]")
      }

      ewd
    }
  }

  private def publishWorkflowSourceToMetadata(id: WorkflowId, workflowSource: WorkflowSource): Unit = {
    val event = MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Workflow), MetadataValue(workflowSource))
    serviceRegistryActor ! PutMetadataAction(event)
  }

  private def pushNamespaceMetadata(validatedNamespace: ValidatedWomNamespace): Unit = {
    val importsMetadata = importedFilesMetadata(validatedNamespace.importedFileContent)
    val wfInputsMetadataEvents = wfInputsMetadata(validatedNamespace.womValueInputs)
    val wfNameMetadataEvent = wfNameMetadata(validatedNamespace.executable.entryPoint.name)
    serviceRegistryActor ! PutMetadataAction(importsMetadata.toVector ++ wfInputsMetadataEvents :+ wfNameMetadataEvent)
  }

  private def pushLanguageToMetadata(languageName: String, languageVersion: String): Unit = {
    val events = List (
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.LanguageName), MetadataValue(languageName)),
      MetadataEvent(
        MetadataKey(workflowId, None, WorkflowMetadataKeys.LanguageVersionName),
        MetadataValue(languageVersion)
      )
    )
    serviceRegistryActor ! PutMetadataAction(events)
  }

  private def wfInputsMetadata(workflowInputs: Map[OutputPort, WomValue]): Iterable[MetadataEvent] = {
    import cromwell.core.CromwellGraphNode.CromwellEnhancedOutputPort

    workflowInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(MetadataKey(workflowId, None, WorkflowMetadataKeys.Inputs)))
      case inputs =>
        inputs flatMap { case (outputPort, womValue) =>
          val inputName = outputPort.fullyQualifiedName
          womValueToMetadataEvents(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Inputs}:$inputName"), womValue)
        }
    }
  }

  private def importedFilesMetadata(imported: Map[String, String]): Iterable[MetadataEvent] = {
    def metadataEventForImportedFile(uri: String, value: String): MetadataEvent = {
      import WorkflowMetadataKeys._
      import cromwell.core.simpleton.WomValueSimpleton._
      // This should only be called on namespaces that are known to have a defined `importUri` so the .get is safe.
      val escapedUri = uri.escapeMeta
      MetadataEvent(
        MetadataKey(workflowId, None, SubmissionSection, SubmissionSection_Imports, escapedUri),
        MetadataValue(value)
      )
    }
    imported map { case (uri, value) => metadataEventForImportedFile(uri, value) }
  }

  private def wfNameMetadata(name: String): MetadataEvent = {
    // Workflow name:
    MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Name), MetadataValue(name))
  }

  private def convertJsonToLabels(json: String): Labels = {
    json.parseJson match {
      case JsObject(inputs) => Labels(inputs.toVector.collect({
        case (key, JsString(value)) => Label(key, value)
      }))
      case _ => Labels(Vector.empty)
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      womNamespace: ValidatedWomNamespace,
                                      workflowOptions: WorkflowOptions,
                                      labels: Labels,
                                      conf: Config,
                                      pathBuilders: List[PathBuilder],
                                      outputRuntimeExtractor: Option[WomOutputRuntimeExtractor]): ErrorOr[EngineWorkflowDescriptor] = {
    val taskCalls = womNamespace.executable.graph.allNodes collect { case taskNode: CommandCallNode => taskNode }
    val defaultBackendName = conf.as[Option[String]]("backend.default")

    val failureModeValidation = validateWorkflowFailureMode(workflowOptions, conf)
    val backendAssignmentsValidation = validateBackendAssignments(taskCalls, workflowOptions, defaultBackendName)

    val callCachingModeValidation = validateCallCachingMode(workflowOptions, conf)

    (failureModeValidation, backendAssignmentsValidation, callCachingModeValidation) mapN {
      case (failureMode, backendAssignments, callCachingMode) =>
        val callable = womNamespace.executable.entryPoint
        val backendDescriptor = BackendWorkflowDescriptor(id, callable, womNamespace.womValueInputs, workflowOptions, labels, hogGroup, List.empty, outputRuntimeExtractor)
        EngineWorkflowDescriptor(callable, backendDescriptor, backendAssignments, failureMode, pathBuilders, callCachingMode)
    }
  }

  private def validateBackendAssignments(calls: Set[CommandCallNode], workflowOptions: WorkflowOptions, defaultBackendName: Option[String]): ErrorOr[Map[CommandCallNode, String]] = {
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
  private def assignBackendUsingRuntimeAttrs(call: CommandCallNode): Option[String] = {
    val runtimeAttributesMap = call.callable.runtimeAttributes.attributes
    runtimeAttributesMap.get(RuntimeBackendKey) map { wdlExpr => evaluateBackendNameExpression(call.fullyQualifiedName, wdlExpr) }
  }

  private def evaluateBackendNameExpression(callName: String, backendNameAsExp: WomExpression): String = {
    backendNameAsExp.evaluateValue(Map.empty, NoIoFunctionSet) match {
      case Valid(runtimeString: WomString) => runtimeString.valueString
      case Valid(x: WomValue) =>
        throw new Exception(s"Non-string values are not currently supported for backends! Cannot use backend '${x.valueString}' to backend to Call: $callName")
      case Invalid(errors) =>
        // TODO WOM: need access to a "source string" for WomExpressions
        // TODO WOM: ErrorOrify this ?
        throw AggregatedMessageException(s"Dynamic backends are not currently supported! Cannot assign backend '$backendNameAsExp' for Call: $callName", errors.toList)
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
