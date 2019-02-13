package cromwell.engine.workflow.lifecycle.materialization

import akka.actor.{ActorRef, FSM, LoggingFSM, Props, Status}
import akka.pattern.pipe
import cats.Monad
import cats.data.EitherT._
import cats.data.Validated.{Invalid, Valid}
import cats.data.NonEmptyList
import cats.effect.IO
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.Checked
import common.exception.{AggregatedMessageException, MessageAggregation}
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr._
import common.validation.Parse._
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowOptions.{ReadFromCache, WorkflowOption, WriteToCache}
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.core.io.AsyncIo
import cromwell.core.labels.{Label, Labels}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.path.PathBuilder
import cromwell.engine._
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.language.CromwellLanguages
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.languages.util.ImportResolver._
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import net.ceedubs.ficus.Ficus._
import spray.json._
import wom.core.{WorkflowSource, WorkflowUrl}
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
            importLocalFilesystem: Boolean, ioActorProxy: ActorRef): Props = {
    Props(new MaterializeWorkflowDescriptorActor(serviceRegistryActor, workflowId, cromwellBackends, importLocalFilesystem, ioActorProxy)).withDispatcher(EngineDispatcher)
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

    val enabled = conf.as[Option[Boolean]]("call-caching.enabled").getOrElse(false)
    val invalidateBadCacheResults = conf.as[Option[Boolean]]("call-caching.invalidate-bad-cache-results").getOrElse(true)
    if (enabled) {
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

      for {
        maybePrefixes <- workflowOptions.getVectorOfStrings("call_cache_hit_path_prefixes")
        callCachingOptions = CallCachingOptions(invalidateBadCacheResults, maybePrefixes)
        mode <- errorOrCallCachingMode(callCachingOptions)
      } yield mode
    }
    else {
      CallCachingOff.validNel
    }
  }
}

// TODO WOM: need to decide where to draw the line between language specific initialization and WOM
class MaterializeWorkflowDescriptorActor(serviceRegistryActor: ActorRef,
                                         workflowId: WorkflowId,
                                         cromwellBackends: => CromwellBackends,
                                         importLocalFilesystem: Boolean,
                                         ioActorProxy: ActorRef) extends LoggingFSM[MaterializeWorkflowDescriptorActorState, Unit] with LazyLogging with WorkflowLogging {

  import MaterializeWorkflowDescriptorActor._
  val tag = self.path.name

  override lazy val workflowIdForLogging = workflowId.toPossiblyNotRoot
  override lazy val rootWorkflowIdForLogging = workflowId.toRoot

  val iOExecutionContext = context.system.dispatchers.lookup("akka.dispatchers.io-dispatcher")
  implicit val ec = context.dispatcher

  startWith(ReadyToMaterializeState, ())

  when(ReadyToMaterializeState) {
    case Event(MaterializeWorkflowDescriptorCommand(workflowSourceFiles, conf), _) =>
      val replyTo = sender()

      workflowOptionsAndPathBuilders(workflowSourceFiles) match {
        case Valid((workflowOptions, pathBuilders)) =>
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
        case Invalid(error) =>
          workflowInitializationFailed(error, replyTo)
          goto(MaterializationFailedState)
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

  private def workflowOptionsAndPathBuilders(sourceFiles: WorkflowSourceFilesCollection): ErrorOr[(WorkflowOptions, Future[List[PathBuilder]])] = {
    val workflowOptionsValidation = validateWorkflowOptions(sourceFiles.workflowOptionsJson)
    workflowOptionsValidation map { workflowOptions =>
      val pathBuilders = EngineFilesystems.pathBuildersForWorkflow(workflowOptions)(context.system)
      (workflowOptions, pathBuilders)
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFilesCollection,
                                      conf: Config,
                                      workflowOptions: WorkflowOptions,
                                      pathBuilders: List[PathBuilder],
                                      engineIoFunctions: EngineIoFunctions): Parse[EngineWorkflowDescriptor] = {

    def findWorkflowSource(workflowSource: Option[WorkflowSource],
                           workflowUrl: Option[WorkflowUrl],
                           resolvers: List[ImportResolver]): Checked[(WorkflowSource, List[ImportResolver])] = {
      (workflowSource, workflowUrl) match {
        case (Some(source), None) => (source, resolvers).validNelCheck
        case (None, Some(url)) =>
          val compoundImportResolver: CheckedAtoB[ImportResolutionRequest, ResolvedImportBundle] = CheckedAtoB.firstSuccess(resolvers.map(_.resolver), s"resolve workflowUrl '$url'")
          val wfSourceAndResolvers: Checked[ResolvedImportBundle] = compoundImportResolver.run(ImportResolutionRequest(url, resolvers))
          wfSourceAndResolvers map { v => (v.source, v.newResolvers) }
        case (Some(_), Some(_)) => "Both workflow source and url can't be supplied".invalidNelCheck
        case (None, None) => "Either workflow source or url has to be supplied".invalidNelCheck
      }
    }

    def findFactory(workflowSource: WorkflowSource): ErrorOr[LanguageFactory] = {
      def chooseFactory(factories: List[LanguageFactory]): Option[LanguageFactory] = factories.find(_.looksParsable(workflowSource))

      val factory: ErrorOr[LanguageFactory] = sourceFiles.workflowType match {
        case Some(languageName) if CromwellLanguages.instance.languages.contains(languageName.toUpperCase) =>
          val language = CromwellLanguages.instance.languages(languageName.toUpperCase)
          sourceFiles.workflowTypeVersion match {
            case Some(v) if language.allVersions.contains(v) => language.allVersions(v).valid
            case Some(other) => s"Unknown version '$other' for workflow language '$languageName'".invalidNel
            case _ => chooseFactory(language.allVersions.values.toList).getOrElse(language.default).valid
          }
        case Some(other) => s"Unknown workflow type: $other".invalidNel[LanguageFactory]
        case None =>
          val allFactories = CromwellLanguages.instance.languages.values.flatMap(_.allVersions.values)
          chooseFactory(allFactories.toList).getOrElse(CromwellLanguages.instance.default.default).validNel
      }

      factory foreach { validFactory =>
        workflowLogger.info(s"Parsing workflow as ${validFactory.languageName} ${validFactory.languageVersionName}")
        pushLanguageToMetadata(validFactory.languageName, validFactory.languageVersionName)
      }
      
      factory
    }

    def buildValidatedNamespace(factory: LanguageFactory, workflowSource: WorkflowSource, importResolvers: List[ImportResolver]): Parse[ValidatedWomNamespace] = {
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

    val zippedResolverCheck: Parse[Option[ImportResolver]] = fromEither[IO](sourceFiles.importsZipFileOption match {
      case None => None.validNelCheck
      case Some(zipContent) => zippedImportResolver(zipContent).toEither.map(Option.apply)
    })

    val labels = convertJsonToLabels(sourceFiles.labelsJson)

    for {
      _ <- publishLabelsToMetadata(id, labels)
      zippedImportResolver <- zippedResolverCheck
      importResolvers = zippedImportResolver.toList ++ localFilesystemResolvers :+ HttpResolver(None, Map.empty)
      sourceAndResolvers <- fromEither[IO](findWorkflowSource(sourceFiles.workflowSource, sourceFiles.workflowUrl, importResolvers))
      _ = if(sourceFiles.workflowUrl.isDefined) publishWorkflowSourceToMetadata(id, sourceAndResolvers._1)
      factory <- errorOrParse(findFactory(sourceAndResolvers._1))
      outputRuntimeExtractor <- errorOrParse(factory.womOutputRuntimeExtractor.toValidated)
      validatedNamespace <- buildValidatedNamespace(factory, sourceAndResolvers._1, sourceAndResolvers._2)
      closeResult = sourceAndResolvers._2.traverse(_.cleanupIfNecessary())
      _ = pushNamespaceMetadata(validatedNamespace)
      ewd <- fromEither[IO](buildWorkflowDescriptor(id, validatedNamespace, workflowOptions, labels, conf, pathBuilders, outputRuntimeExtractor).toEither)
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

  private def publishLabelsToMetadata(rootWorkflowId: WorkflowId, labels: Labels): Parse[Unit] = {
    val defaultLabel = "cromwell-workflow-id" -> s"cromwell-$rootWorkflowId"
    val customLabels = labels.asMap
    Monad[Parse].pure(labelsToMetadata(customLabels + defaultLabel, rootWorkflowId))
  }

  protected def labelsToMetadata(labels: Map[String, String], workflowId: WorkflowId): Unit = {
    labels foreach { case (k, v) =>
      serviceRegistryActor ! PutMetadataAction(MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Labels}:$k"), MetadataValue(v)))
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
        val backendDescriptor = BackendWorkflowDescriptor(id, callable, womNamespace.womValueInputs, workflowOptions, labels, outputRuntimeExtractor)
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

  private def validateWorkflowOptions(workflowOptions: String): ErrorOr[WorkflowOptions] = {
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
