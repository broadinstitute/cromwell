package cromwell.engine.workflow

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.{Dispatcher, WorkflowId}
import cromwell.database.sql.SqlDatabase
import cromwell.database.sql.tables.DockerHashStoreEntry
import cromwell.docker.DockerHashActor.{DockerHashFailureResponse, DockerHashSuccessResponse}
import cromwell.docker.{DockerClientHelper, DockerHashRequest, DockerHashResult, DockerImageIdentifier}
import cromwell.engine.workflow.WorkflowActor.{RestartExistingWorkflow, StartMode}
import cromwell.engine.workflow.WorkflowDockerLookupActor._
import cromwell.services.SingletonServicesStore
import cromwell.core.Dispatcher.EngineDispatcher

import lenthall.util.TryUtil

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

/**
  * Ensures docker hash consistency throughout a workflow.
  *
  * Caches successful docker hash lookups and serves them to subsequent identical requests.
  * Persists those hashes in the database to be resilient to server restarts.
  *
  * Failure modes:
  * 1) Failure to load hashes from the DB upon restart.
  * 2) Failure to parse hashes from the DB upon restart.
  * 3) Failure to write a hash result to the DB.
  * 4) Failure to look up a docker hash.
  * 5) Timeout from DockerHashActor.
  *
  * Behavior:
  * 1-2) Return a terminal lookup failure for all pending requests, transition to a permanently Failed state in which any
  *      future requests will immediately return lookup failure.  The JobPreparation actor should fail in response to this
  *      lookup termination message, which in turn should fail the workflow.
  * 3-5) Return lookup failure for the current request and all pending requests for the same tag.  Any future requests
  *      for this tag will be attempted again.
  */

class WorkflowDockerLookupActor private[workflow](workflowId: WorkflowId, val dockerHashingActor: ActorRef, startMode: StartMode, databaseInterface: SqlDatabase)
  extends LoggingFSM[WorkflowDockerLookupActorState, WorkflowDockerLookupActorData] with DockerClientHelper {

  implicit val ec = context.system.dispatchers.lookup(Dispatcher.EngineDispatcher)

  // Amount of time after which the docker request should be considered lost and sent again.
  override protected def backpressureTimeout: FiniteDuration = 10 seconds
  // A multiplier for the amount of time to wait when we get a Backpressure response before sending the request again.
  // This effectively bounds the jitter.
  override protected def backpressureRandomizerFactor: Double = 0.5D

  context.become(dockerReceive orElse receive)

  startWith(
    stateName = if (startMode == RestartExistingWorkflow) AwaitingFirstRequestOnRestart else Running,
    stateData = WorkflowDockerLookupActorData.empty
  )

  // `AwaitingFirstRequestOnRestart` is only used in restart scenarios.  This state waits until there's at least one hash
  // request before trying to load the docker hash mappings.  This is so we'll have at least one `JobPreparationActor`
  // reference available to message with a terminal failure in case the reading or parsing of these mappings fails.
  when(AwaitingFirstRequestOnRestart) {
    case Event(request: DockerHashRequest, data) =>
      loadDockerHashStoreEntries()
      goto(LoadingCache) using data.addHashRequest(request, sender())
  }

  // Waiting for a response from the database with the hash mapping for this workflow.
  when(LoadingCache) {
    case Event(DockerHashStoreLoadingSuccess(dockerHashEntries), data) =>
      loadCacheAndHandleHashRequests(dockerHashEntries, data)
    case Event(request: DockerHashRequest, data) =>
      stay using data.addHashRequest(request, sender())
  }

  // This is the normal operational mode.
  when(Running) {
    // This tag has already been looked up and its hash is in the mappings cache.
    case Event(request: DockerHashRequest, data) if data.mappings.contains(request.dockerImageID) =>
      sender ! DockerHashSuccessResponse(data.mappings(request.dockerImageID), request)
      stay()
    // A request for the hash for this tag has already been made to the hashing actor.  Don't request the hash again,
    // just add this sender to the list of replyTos for when the hash arrives.
    case Event(request: DockerHashRequest, data) if data.hashRequests.contains(request) =>
      stay using data.addHashRequest(request, sender())
    // This tag has not (successfully) been looked up before, so look it up now.
    case Event(request: DockerHashRequest, data) =>
      requestDockerHash(request, data)
    case Event(dockerResponse: DockerHashSuccessResponse, data) =>
      persistDockerHash(dockerResponse, data)
      stay()
    case Event(dockerResponse: DockerHashFailureResponse, data) =>
      handleLookupFailure(dockerResponse, data)
    case Event(DockerHashStoreSuccess(response), data) =>
      recordMappingAndRespond(response, data)
    case Event(DockerHashStoreFailure(request, e), data) =>
      handleStoreFailure(request, new Exception(s"Failure storing docker hash for ${request.dockerImageID.fullName}", e), data)
  }

  // In state Terminal we reject all requests with the cause set in the state data.
  when(Terminal) {
    case Event(request: DockerHashRequest, data) =>
      sender() ! WorkflowDockerLookupFailure(data.failureCause.orNull, request)
      stay()
  }

  private def fail(reason: Throwable): State = {
    self ! TransitionToFailed(reason)
    stay()
  }

  whenUnhandled {
    case Event(DockerHashActorTimeout(request), data) =>
      val reason = new Exception(s"Timeout looking up docker hash")
      data.hashRequests(request) foreach { _ ! WorkflowDockerLookupFailure(reason, request) }
      val updatedData = data.copy(hashRequests = data.hashRequests - request)
      stay() using updatedData
    case Event(TransitionToFailed(cause), data) =>
      log.error(cause, s"Workflow Docker lookup actor for $workflowId transitioning to Failed")
      val updatedData = respondToAllRequestsWithTerminalFailure(FailedException, data)
      goto(Terminal) using updatedData.withFailureCause(FailedException)
  }

  /**
    * Load mappings from the database into the state data, reply to queued requests which have mappings, and initiate
    * hash lookups for requests which don't have mappings.
    */
  private def loadCacheAndHandleHashRequests(hashEntries: Map[String, String], data: WorkflowDockerLookupActorData): State = {
    val dockerMappingsTry = hashEntries map {
      case (dockerTag, dockerHash) => DockerImageIdentifier.fromString(dockerTag) -> DockerHashResult.fromString(dockerHash)
    }

    TryUtil.sequenceKeyValues(dockerMappingsTry) match {
      case Success(dockerMappings) =>
        // Figure out which of the queued requests already have established mappings.
        val (hasMappings, doesNotHaveMappings) = data.hashRequests.partition { case (request, _) => dockerMappings.contains(request.dockerImageID) }

        // The requests which have mappings receive success responses.
        hasMappings foreach { case (request, replyTos) =>
          val result = dockerMappings(request.dockerImageID)
          replyTos foreach { _ ! DockerHashSuccessResponse(result, request)}
        }

        // The requests without mappings need to be looked up.
        doesNotHaveMappings.keys foreach { sendDockerCommand(_) }

        // Update state data accordingly.
        val newData = data.copy(hashRequests = doesNotHaveMappings, mappings = dockerMappings, failureCause = None)
        goto(Running) using newData

      case Failure(e) =>
        fail(new Exception("Failed to parse docker tag -> hash mappings from DB", e))
    }
  }

  private def requestDockerHash(request: DockerHashRequest, data: WorkflowDockerLookupActorData): State = {
    sendDockerCommand(request)
    val replyTo = sender()
    val updatedData = data.copy(hashRequests = data.hashRequests + (request -> List(replyTo)))
    stay using updatedData
  }

  private def recordMappingAndRespond(response: DockerHashSuccessResponse, data: WorkflowDockerLookupActorData): State = {
    // Add the new label to hash mapping to the current set of mappings.
    val request = response.request
    data.hashRequests(request) foreach { _ ! DockerHashSuccessResponse(response.dockerHash, request) }
    val updatedData = data.copy(hashRequests = data.hashRequests - request, mappings = data.mappings + (request.dockerImageID -> response.dockerHash))
    stay using updatedData
  }

  private def respondToAllRequests(reason: Throwable,
                                   data: WorkflowDockerLookupActorData,
                                   messageBuilder: (Throwable, DockerHashRequest) => WorkflowDockerLookupResponse): WorkflowDockerLookupActorData = {
    data.hashRequests foreach { case (request, replyTos) =>
      replyTos foreach { _ ! messageBuilder(reason, request) }
    }
    data.clearHashRequests
  }

  private def respondToAllRequestsWithTerminalFailure(reason: Throwable, data: WorkflowDockerLookupActorData): WorkflowDockerLookupActorData = {
    respondToAllRequests(reason, data, WorkflowDockerTerminalFailure.apply)
  }

  private def persistDockerHash(response: DockerHashSuccessResponse, data: WorkflowDockerLookupActorData): Unit = {
    val dockerHashStoreEntry = DockerHashStoreEntry(workflowId.toString, response.request.dockerImageID.fullName, response.dockerHash.algorithmAndHash)
    databaseInterface.addDockerHashStoreEntry(dockerHashStoreEntry) onComplete {
      case Success(_) => self ! DockerHashStoreSuccess(response)
      case Failure(ex) => self ! DockerHashStoreFailure(response.request, ex)
    }
  }

  private def handleLookupFailure(dockerResponse: DockerHashFailureResponse, data: WorkflowDockerLookupActorData): State = {
    // Fail all pending requests.  This logic does not blacklist the tag, which will allow lookups to be attempted
    // again in the future.
    val failureResponse = WorkflowDockerLookupFailure(new Exception(dockerResponse.reason), dockerResponse.request)
    val request = dockerResponse.request
    data.hashRequests(request) foreach { _ ! failureResponse }

    val updatedData = data.copy(hashRequests = data.hashRequests - request)
    stay using updatedData
  }

  private def handleStoreFailure(dockerHashRequest: DockerHashRequest, reason: Throwable, data: WorkflowDockerLookupActorData): State = {
    data.hashRequests(dockerHashRequest) foreach { _ ! WorkflowDockerLookupFailure(reason, dockerHashRequest) }
    // Remove these requesters from the collection of those awaiting hashes.
    stay() using data.copy(hashRequests = data.hashRequests - dockerHashRequest)
  }

  def loadDockerHashStoreEntries(): Unit = {
    databaseInterface.queryDockerHashStoreEntries(workflowId.toString) onComplete {
      case Success(dockerHashEntries) =>
        val dockerMappings = dockerHashEntries.map(entry => entry.dockerTag -> entry.dockerHash).toMap
        self ! DockerHashStoreLoadingSuccess(dockerMappings)
      case Failure(ex) =>
        fail(new RuntimeException("Failed to load docker tag -> hash mappings from DB", ex))
    }
  }

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    message match {
      case r: DockerHashRequest => self ! DockerHashActorTimeout(r)
    }
  }
}

object WorkflowDockerLookupActor {
  /* States */
  sealed trait WorkflowDockerLookupActorState
  case object AwaitingFirstRequestOnRestart extends WorkflowDockerLookupActorState
  case object LoadingCache extends WorkflowDockerLookupActorState
  case object Running extends WorkflowDockerLookupActorState
  case object Terminal extends WorkflowDockerLookupActorState
  private val FailedException =
    new Exception(s"The WorkflowDockerLookupActor has failed. Subsequent docker tags for this workflow will not be resolved.")

  /* Internal ADTs */
  final case class DockerRequestContext(dockerHashRequest: DockerHashRequest, replyTo: ActorRef)
  sealed trait DockerHashStoreResponse
  final case class DockerHashStoreSuccess(successResponse: DockerHashSuccessResponse) extends DockerHashStoreResponse
  final case class DockerHashStoreFailure(dockerHashRequest: DockerHashRequest, reason: Throwable) extends DockerHashStoreResponse
  final case class DockerHashStoreLoadingSuccess(dockerMappings: Map[String, String])
  final case class DockerHashActorTimeout(request: DockerHashRequest)

  /* Messages */
  sealed trait WorkflowDockerLookupActorMessage
  private final case class TransitionToFailed(cause: Throwable) extends WorkflowDockerLookupActorMessage

  /* Responses */
  sealed trait WorkflowDockerLookupResponse
  final case class WorkflowDockerLookupFailure(reason: Throwable, request: DockerHashRequest) extends WorkflowDockerLookupResponse
  final case class WorkflowDockerTerminalFailure(reason: Throwable, request: DockerHashRequest) extends WorkflowDockerLookupResponse

  def props(workflowId: WorkflowId, dockerHashingActor: ActorRef, startMode: StartMode, databaseInterface: SqlDatabase = SingletonServicesStore.databaseInterface) = {
    Props(new WorkflowDockerLookupActor(workflowId, dockerHashingActor, startMode, databaseInterface)).withDispatcher(EngineDispatcher)
  }

  object WorkflowDockerLookupActorData {
    def empty = WorkflowDockerLookupActorData(hashRequests = Map.empty, mappings = Map.empty, failureCause = None)
  }

  final case class WorkflowDockerLookupActorData(hashRequests: Map[DockerHashRequest, List[ActorRef]],
                                                 mappings: Map[DockerImageIdentifier, DockerHashResult],
                                                 failureCause: Option[Throwable]) {
    /**
      * Add the specified request and replyTo to this state data.
      *
      * @param request The request to be added.
      * @param replyTo The actor to be informed of the hash or the failure to look up the hash.
      * @return State data with the added request and replyTo.
      */
    def addHashRequest(request: DockerHashRequest, replyTo: ActorRef): WorkflowDockerLookupActorData = {
      // Prepend this `ActorRef` to the list of `ActorRef`s awaiting the hash for this request, or to Nil if this is the first.
      val alreadyAwaiting = hashRequests.getOrElse(request, Nil)
      this.copy(hashRequests = hashRequests + (request -> (replyTo :: alreadyAwaiting)))
    }

    /**
      * Empty the collection of hash requests.
      * @return State data with all awaiting hash requests removed.
      */
    def clearHashRequests: WorkflowDockerLookupActorData = this.copy(hashRequests = Map.empty)

    /**
      * Add this failure cause to the state data.
      * @param cause The failure cause.
      * @return Updated state data.
      */
    def withFailureCause(cause: Throwable): WorkflowDockerLookupActorData = this.copy(failureCause = Option(cause))
  }
}
