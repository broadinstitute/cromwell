package cromwell.engine.workflow.lifecycle.materialization

import akka.actor.{ActorRef, FSM, LoggingFSM, Props, Status}
import akka.pattern.pipe
import better.files.File
import cats.Monad
import cats.data.EitherT._
import cats.data.Validated.{Invalid, Valid}
import cats.data.{EitherT, NonEmptyList}
import cats.effect.IO
import cats.instances.vector._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.Checked
import common.exception.{AggregatedMessageException, MessageAggregation}
import common.validation.Checked._
import common.validation.ErrorOr._
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.CromwellGraphNode._
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
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cwl.CwlDecoder
import cwl.CwlDecoder.Parse
import net.ceedubs.ficus.Ficus._
import spray.json._
import wdl._
import wom.core.WorkflowSource
import wom.executable.Executable
import wom.executable.Executable.ResolvedExecutableInputs
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{Graph, TaskCallNode}
import wom.values.{WomSingleFile, WomString, WomValue}

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
  case object MaterializingState extends MaterializeWorkflowDescriptorActorState
  case object MaterializationSuccessfulState extends MaterializeWorkflowDescriptorActorTerminalState
  case object MaterializationFailedState extends MaterializeWorkflowDescriptorActorTerminalState
  case object MaterializationAbortedState extends MaterializeWorkflowDescriptorActorTerminalState

  private val DefaultWorkflowFailureMode = NoNewCalls.toString

  /*
    * Internal ADT
   */
  private case class ValidatedWomNamespace(executable: Executable, graph: Graph, womValueInputs: Map[OutputPort, WomValue])

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

      (readFromCache, writeToCache) mapN {
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

// TODO WOM: need to decide where to draw the line between language specific initialization and WOM
class MaterializeWorkflowDescriptorActor(serviceRegistryActor: ActorRef,
                                         val workflowIdForLogging: WorkflowId,
                                         cromwellBackends: => CromwellBackends,
                                         importLocalFilesystem: Boolean) extends LoggingFSM[MaterializeWorkflowDescriptorActorState, Unit] with LazyLogging with WorkflowLogging {

  import MaterializeWorkflowDescriptorActor._

  val tag = self.path.name

  val iOExecutionContext = context.system.dispatchers.lookup("akka.dispatchers.io-dispatcher")
  implicit val ec = context.dispatcher

  startWith(ReadyToMaterializeState, ())

  when(ReadyToMaterializeState) {
    case Event(MaterializeWorkflowDescriptorCommand(workflowSourceFiles, conf), _) =>
      val replyTo = sender()

      workflowOptionsAndPathBuilders(workflowSourceFiles) match {
        case Valid((workflowOptions, pathBuilders)) =>
          val futureDescriptor: Future[ErrorOr[EngineWorkflowDescriptor]] = pathBuilders flatMap {
            buildWorkflowDescriptor(workflowIdForLogging, workflowSourceFiles, conf, workflowOptions, _).
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
      val pathBuilders = EngineFilesystems.pathBuildersForWorkflow(workflowOptions)(context.system, context.dispatcher)
      (workflowOptions, pathBuilders)
    }
  }

  private def buildWorkflowDescriptor(id: WorkflowId,
                                      sourceFiles: WorkflowSourceFilesCollection,
                                      conf: Config,
                                      workflowOptions: WorkflowOptions,
                                      pathBuilders: List[PathBuilder]): Parse[EngineWorkflowDescriptor] = {
    val namespaceValidation: Parse[ValidatedWomNamespace] = sourceFiles.workflowType match {
      case Some(wdl) if wdl.equalsIgnoreCase("wdl") => fromEither[IO](validateWdlNamespace(sourceFiles, workflowOptions, pathBuilders).toEither)
      case Some(cwl) if cwl.equalsIgnoreCase("cwl") => validateCwlNamespace(sourceFiles, workflowOptions, pathBuilders)
      case Some(other) => fromEither[IO](s"Unknown workflow type: $other".invalidNelCheck[ValidatedWomNamespace])
      case None => fromEither[IO]("Need a workflow type here !".invalidNelCheck[ValidatedWomNamespace])
    }
    val labelsValidation: Parse[Labels] = fromEither[IO](validateLabels(sourceFiles.labelsJson).toEither)

    for {
      validatedNamespace <- namespaceValidation
      labels <- labelsValidation
      _ <- pushWfNameMetadataService(validatedNamespace.executable.entryPoint.name)
      _ <- publishLabelsToMetadata(id, validatedNamespace.executable.entryPoint.name, labels)
      ewd <- fromEither[IO](buildWorkflowDescriptor(id, sourceFiles, validatedNamespace, workflowOptions, labels, conf, pathBuilders).toEither)
    } yield ewd
  }


  private def pushWfNameMetadataService(name: String): Parse[Unit] = {
    // Workflow name:
    val nameEvent = MetadataEvent(MetadataKey(workflowIdForLogging, None, WorkflowMetadataKeys.Name), MetadataValue(name))

    Monad[Parse].pure(serviceRegistryActor ! PutMetadataAction(nameEvent))
  }

  private def publishLabelsToMetadata(rootWorkflowId: WorkflowId, unqualifiedName: String, labels: Labels): Parse[Unit] = {
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
                                      sourceFiles: WorkflowSourceFilesCollection,
                                      womNamespace: ValidatedWomNamespace,
                                      workflowOptions: WorkflowOptions,
                                      labels: Labels,
                                      conf: Config,
                                      pathBuilders: List[PathBuilder]): ErrorOr[EngineWorkflowDescriptor] = {
    val taskCalls = womNamespace.graph.allNodes collect { case taskNode: TaskCallNode => taskNode }
    val defaultBackendName = conf.as[Option[String]]("backend.default")

    val failureModeValidation = validateWorkflowFailureMode(workflowOptions, conf)
    val backendAssignmentsValidation = validateBackendAssignments(taskCalls, workflowOptions, defaultBackendName)

    val callCachingModeValidation = validateCallCachingMode(workflowOptions, conf)

    (failureModeValidation, backendAssignmentsValidation, callCachingModeValidation) mapN {
      case (failureMode, backendAssignments, callCachingMode) =>
        val callable = womNamespace.executable.entryPoint
        val backendDescriptor = BackendWorkflowDescriptor(id, callable, womNamespace.womValueInputs, workflowOptions, labels)
        EngineWorkflowDescriptor(callable, backendDescriptor, backendAssignments, failureMode, pathBuilders, callCachingMode)
    }
  }

  private def pushWfInputsToMetadataService(workflowInputs: Map[OutputPort, WomValue]): Unit = {
    // Inputs
    val inputEvents = workflowInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(MetadataKey(workflowIdForLogging, None,WorkflowMetadataKeys.Inputs)))
      case inputs =>
        inputs flatMap { case (outputPort, womValue) =>
          val inputName = outputPort.fullyQualifiedName
          womValueToMetadataEvents(MetadataKey(workflowIdForLogging, None, s"${WorkflowMetadataKeys.Inputs}:$inputName"), womValue)
        }
    }

    serviceRegistryActor ! PutMetadataAction(inputEvents)
  }

  private def validateBackendAssignments(calls: Set[TaskCallNode], workflowOptions: WorkflowOptions, defaultBackendName: Option[String]): ErrorOr[Map[TaskCallNode, String]] = {
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
  private def assignBackendUsingRuntimeAttrs(call: TaskCallNode): Option[String] = {
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

  // TODO WOM: resurect ?
  //  private def validateDeclarations(namespace: WdlNamespaceWithWorkflow,
  //                                   options: WorkflowOptions,
  //                                   coercedInputs: WorkflowCoercedInputs,
  //                                   pathBuilders: List[PathBuilder]): ErrorOr[WorkflowCoercedInputs] = {
  //    namespace.staticDeclarationsRecursive(coercedInputs, NoFunctions) match {
  //      case Success(d) => d.validNel
  //      case Failure(e) => s"Workflow has invalid declarations: ${e.getMessage}".invalidNel
  //    }
  //  }

  private def validateImportsDirectory(zipContents: Array[Byte], parentPath: Option[Path] = None): ErrorOr[Path] = {

    def makeZipFile: Try[Path] = Try {
      DefaultPathBuilder.createTempFile("", ".zip", parentPath).writeByteArray(zipContents)(OpenOptions.default)
    }

    def unZipFile(f: Path) = Try {
      val unzippedFile = f.unzipTo(parentPath)
      val unzippedFileContents = unzippedFile.list.toSeq.head
      if (unzippedFileContents.isDirectory) unzippedFileContents else unzippedFile
    }

    val importsFile = for {
      zipFile <- makeZipFile
      unzipped <- unZipFile(zipFile)
      _ <- Try(zipFile.delete(swallowIOExceptions = true))
    } yield unzipped

    importsFile match {
      case Success(unzippedDirectory: Path) => unzippedDirectory.validNel
      case Failure(t) => t.getMessage.invalidNel
    }
  }

  private def validateWorkflowNameLengths(namespace: WdlNamespaceWithWorkflow): Checked[Unit] = {
    import common.validation.Checked._
    def allWorkflowNames(n: WdlNamespace): Seq[String] = n.workflows.map(_.unqualifiedName) ++ n.namespaces.flatMap(allWorkflowNames)
    val tooLong = allWorkflowNames(namespace).filter(_.length >= 100)
    if (tooLong.nonEmpty) {
      ("Workflow names must be shorter than 100 characters: " + tooLong.mkString(" ")).invalidNelCheck
    } else {
      ().validNelCheck
    }
  }

  private def validateCwlNamespace(source: WorkflowSourceFilesCollection,
                                   workflowOptions: WorkflowOptions,
                                   pathBuilders: List[PathBuilder]): Parse[ValidatedWomNamespace] = {
    // TODO WOM: CwlDecoder takes a file so write it to disk for now

    val cwlFile: File = File.newTemporaryFile(prefix = workflowIdForLogging.toString).write(source.workflowSource)

    def unzipDependencies: Parse[Unit]  = source match {
      case wsfwdz: WorkflowSourceFilesWithDependenciesZip =>
        EitherT.
          fromEither[IO]{validateImportsDirectory(wsfwdz.importsZip, Some(DefaultPathBuilder.build(cwlFile.parent.path))).void.toEither}
      case _ => Monad[Parse].unit
    }

    import cwl.AcceptAllRequirements
    for {
      _ <- unzipDependencies
      cwl <- CwlDecoder.decodeAllCwl(cwlFile, source.workflowRoot)
      executable <- fromEither[IO](cwl.womExecutable(AcceptAllRequirements, Option(source.inputsJson)))
      ioFunctions = new WdlFunctions(pathBuilders)
      validatedWomNamespace <- fromEither[IO](validateWomNamespace(executable, ioFunctions))
    } yield validatedWomNamespace
  }

  private def validateWdlNamespace(source: WorkflowSourceFilesCollection,
                                   workflowOptions: WorkflowOptions,
                                   pathBuilders: List[PathBuilder]): ErrorOr[ValidatedWomNamespace] = {
    import cats.instances.either._
    import cats.instances.list._
    import cats.syntax.either._
    import cats.syntax.functor._
    import common.validation.Checked._

    def checkTypes(namespace: WdlNamespaceWithWorkflow, inputs: Map[OutputPort, WomValue]): Checked[Unit] = {
      val allDeclarations = namespace.workflow.declarations ++ namespace.workflow.calls.flatMap(_.declarations)
      val list: List[Checked[Unit]] = inputs.map({ case (k, v) =>
        allDeclarations.find(_.fullyQualifiedName == k) match {
          case Some(decl) if decl.womType.coerceRawValue(v).isFailure =>
            s"Invalid right-side type of '$k'.  Expecting ${decl.womType.toDisplayString}, got ${v.womType.toDisplayString}".invalidNelCheck[Unit]
          case _ => ().validNelCheck
        }
      }).toList

      list.sequence[Checked, Unit].void
    }

    val baseResolvers: List[String => WorkflowSource] = if (importLocalFilesystem) {
      List(WdlNamespace.fileResolver, WdlNamespace.httpResolver)
    } else {
      List(WdlNamespace.httpResolver)
    }

    import common.validation.Validation._

    val wdlNamespaceValidation = source match {
      case w: WorkflowSourceFilesWithDependenciesZip =>
        for {
          importsDir <- validateImportsDirectory(w.importsZip)
          betterFilesImportsDir = better.files.File(importsDir.pathAsString)
          directoryResolver = WdlNamespace.directoryResolver(betterFilesImportsDir): String => WorkflowSource
          resolvers = directoryResolver +: baseResolvers
          wf <- WdlNamespaceWithWorkflow.load(w.workflowSource, resolvers).toErrorOr
          _ = importsDir.delete(swallowIOExceptions = true)
        } yield wf
      case w: WorkflowSourceFilesWithoutImports =>
        WdlNamespaceWithWorkflow.load(w.workflowSource, baseResolvers).toErrorOr
    }

    /* Publish `MetadataEvent`s for all imports in this `WdlNamespace`. */
    def publishImportMetadata(wdlNamespace: WdlNamespace): Unit = {
      def metadataEventForImportedNamespace(ns: WdlNamespace): MetadataEvent = {
        import MetadataKey._
        import WorkflowMetadataKeys._
        // This should only be called on namespaces that are known to have a defined `importUri` so the .get is safe.
        val escapedUri = ns.importUri.get.escapeMeta
        MetadataEvent(MetadataKey(
          workflowIdForLogging, None, SubmissionSection, SubmissionSection_Imports, escapedUri), MetadataValue(ns.sourceString))
      }

      // Descend the namespace looking for imports and construct `MetadataEvent`s for them.
      def collectImportEvents: List[MetadataEvent] = {
        wdlNamespace.allNamespacesRecursively flatMap { ns =>
          ns.importUri.toList collect {
            // Do not publish import events for URIs which correspond to literal strings as these are the top-level
            // submitted workflow.
            case uri if uri != WdlNamespace.WorkflowResourceString => metadataEventForImportedNamespace(ns)
          }
        }
      }

      serviceRegistryActor ! PutMetadataAction(collectImportEvents)
    }

    (for {
      wdlNamespace <- wdlNamespaceValidation.toEither
      _ <- validateWorkflowNameLengths(wdlNamespace)
      _ = publishImportMetadata(wdlNamespace)
      womExecutable <- wdlNamespace.womExecutable(Option(source.inputsJson))
      ioFunctions = new WdlFunctions(pathBuilders)
      validatedWomNamespace <- validateWomNamespace(womExecutable, ioFunctions)
      _ <- checkTypes(wdlNamespace, validatedWomNamespace.womValueInputs)
      _ = pushWfInputsToMetadataService(validatedWomNamespace.womValueInputs)
    } yield validatedWomNamespace).toValidated
  }

  private def validateWomNamespace(womExecutable: Executable, ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace] = for {
    evaluatedInputs <- validateExecutableInputs(womExecutable.resolvedExecutableInputs, ioFunctions).toEither
    validatedWomNamespace = ValidatedWomNamespace(womExecutable, womExecutable.graph, evaluatedInputs)
    _ <- validateWdlFiles(validatedWomNamespace.womValueInputs)
  } yield validatedWomNamespace

  /*
    * At this point input values are either a WomValue (if it was provided through the input file)
    * or a WomExpression (if we fell back to the default).
    * We assume that default expressions do NOT reference any "variables" (other inputs, call outputs ...)
    * Should this assumption prove not sufficient InstantiatedExpressions or ExpressionNodes would need to be provided
    * instead so that they can be evaluated JIT.
    * Note that the ioFunctions use engine level pathBuilders. This means that their credentials come from the engine section
    * of the configuration, and are not backend specific.
   */
  private def validateExecutableInputs(inputs: ResolvedExecutableInputs, ioFunctions: IoFunctionSet): ErrorOr[Map[OutputPort, WomValue]] = {
    inputs.traverse[OutputPort, WomValue] {
      case (key, value) => value.fold(ResolvedExecutableInputsPoly).apply(ioFunctions) map { key -> _ }
    }
  }

  private def validateLabels(json: WorkflowJson): ErrorOr[Labels] = {

    def toLabels(inputs: Map[String, JsValue]): ErrorOr[Labels] = {
      val vectorOfValidatedLabel: Vector[ErrorOr[Label]] = inputs.toVector map {
        case (key, JsString(s)) => Label.validateLabel(key, s)
        case (key, other) => s"Invalid label $key: $other : Labels must be strings. ${Label.LabelExpectationsMessage}".invalidNel
      }

      vectorOfValidatedLabel.sequence[ErrorOr, Label] map { validatedVectorofLabel => Labels(validatedVectorofLabel) }
    }

    Try(json.parseJson) match {
      case Success(JsObject(inputs)) => toLabels(inputs)
      case Failure(reason: Throwable) => s"Workflow contains invalid labels JSON: ${reason.getMessage}".invalidNel
      case _ => """Invalid workflow labels JSON. Expected a JsObject of "labelKey": "labelValue" values.""".invalidNel
    }
  }

  private def validateWdlFiles(workflowInputs: Map[OutputPort, WomValue]): Checked[Unit] = {
    
    def prefix(port: OutputPort) = s"Invalid value for File input '${port.fullyQualifiedName}':"
    
    val failedFiles = workflowInputs collect {
      case (port , WomSingleFile(value)) if value.startsWith("\"gs://") =>
        s"""${prefix(port)} $value starts with a '"'"""
      case (port , WomSingleFile(value)) if value.isEmpty =>
        s"${prefix(port)} empty value"
    }

    NonEmptyList.fromList(failedFiles.toList) match {
      case Some(errors) => Left(errors)
      case None => Right(())
    }
  }

  private def validateWorkflowOptions(workflowOptions: WorkflowJson): ErrorOr[WorkflowOptions] = {
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
