package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.BackendJobDescriptor
import cromwell.core.JobKey
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{HashingFailedMessage, PlaceholderHashKeyExpansion, SuccessfulHashResultMessage}
import cromwell.engine.workflow.lifecycle.execution.callcaching.BackendSpecificHasherActor._
import org.apache.commons.lang3.NotImplementedException
import wdl4s.values.WdlFile

/**
  * Used if a backend has not provided a backend hasher actor.
  *
  * Backend specific hasher actors should extend and reply to commands in the same way.
  */
class BackendSpecificHasherActor(mode: CallCachingMode) extends Actor with ActorLogging {

  override def receive: Receive = {
    // We can't do file hashes, so respond with a bunch of failures.
    case JobFileHashRequests(jobKey, fileHashRequests) =>
      fileHashRequests map { fileHashRequest => HashingFailedMessage(fileHashRequest.hashKey, new NotImplementedException("This type of hashing is not implemented for this backend")) } foreach { sender ! _ }

    // Expand the placeholder key, and then send back a failure message for each of them!
    case RuntimeAttributesHashesRequest(placeholderHashKey, jobDescriptor) =>
      val newKeys = Map (
        mode.lookupDockerHashes -> HashKey("runtime attribute: docker(Hash lookup)"),
        mode.hashDockerNames -> HashKey("runtime attribute: docker(Name only)")
      ) collect { case (true, x) => x }
      sender ! PlaceholderHashKeyExpansion(placeholderHashKey, newKeys)
      newKeys foreach { key => sender ! HashingFailedMessage(key, new NotImplementedException("This type of hashing is not implemented for this backend")) }

    // Since all these methods are sync, we'll never have any work outstanding. So we don't need to do anything.
    case cancellation: HashesNoLongerRequired => // Nothing
  }
}

object BackendSpecificHasherActor {

  def props(mode: CallCachingMode): Props = Props(new BackendSpecificHasherActor(mode))

  case class SingleFileHashRequest(hashKey: HashKey, file: WdlFile)

  sealed trait BackendSpecificHasherCommand
  case class JobFileHashRequests(jobKey: JobKey, files: Iterable[SingleFileHashRequest]) extends BackendSpecificHasherCommand
  case class RuntimeAttributesHashesRequest(placeholderHashKey: HashKey, jobDescriptor: BackendJobDescriptor)
  case class HashesNoLongerRequired(jobKey: JobKey)

  sealed trait BackendSpecificHasherResponse extends SuccessfulHashResultMessage
  case class FileHashesResponse(hashResult: HashResult) extends BackendSpecificHasherResponse { override def hashes = List(hashResult) }
}
