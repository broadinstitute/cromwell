package cromwell.backend.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor}
import cromwell.backend.callcaching.BackendSpecificHasherActor._
import cromwell.core.JobKey
import cromwell.core.callcaching._
import wdl4s.values.WdlFile
import cromwell.core.callcaching.HashValue.StringMd5er

/**
  * Handles requests between the EngineJobHasherActor and the BackendHashingMethods.
  */
class BackendSpecificHasherActor(backendHashingMethods: BackendHashingMethods) extends Actor with ActorLogging {

  override def receive: Receive = {
    // We can't do file hashes, so respond with a bunch of failures.
    case JobFileHashRequests(jobKey, fileHashRequests) => fileHashRequests foreach  { backendHashingMethods.fileContentsHasherActor ! _ }

    case s: HashResultMessage => context.parent ! s

    // Expand the placeholder key, and then send back a failure message for each of them!
    case RuntimeAttributesHashesRequest(jobDescriptor) =>
      def toHashKey(attributeName: String) = HashKey("runtime attribute: " + attributeName)
      val hashesNeeded = backendHashingMethods.hashableRuntimeAttributes
      sender ! RuntimeAttributesHashKeyPlaceholderExpansion(hashesNeeded map toHashKey)

      val hashValues = hashesNeeded map { attributeName =>
        HashResult(toHashKey(attributeName), jobDescriptor.runtimeAttributes.get(attributeName).map(_.valueString).map(_.md5HashValue).getOrElse(UnspecifiedRuntimeAttributeHashValue))
      }

      sender ! RuntimeAttributeHashesResponse(hashValues.toSet)

    // Since all these methods are sync, we'll never have any work outstanding. So we don't need to do anything.
    case cancellation: HashesNoLongerRequired =>
      context.stop(self)
  }
}

object BackendSpecificHasherActor {

  def props(backendHashingMethods: BackendHashingMethods): Props = Props(new BackendSpecificHasherActor(backendHashingMethods))

  case class SingleFileHashRequest(hashKey: HashKey, file: WdlFile, initializationData: Option[BackendInitializationData])

  sealed trait BackendSpecificHasherCommand
  case class JobFileHashRequests(jobKey: JobKey, files: Iterable[SingleFileHashRequest]) extends BackendSpecificHasherCommand
  case class RuntimeAttributesHashesRequest(jobDescriptor: BackendJobDescriptor) extends BackendSpecificHasherCommand
  case class HashesNoLongerRequired(jobKey: JobKey)

  sealed trait BackendSpecificHasherResponse extends SuccessfulHashResultMessage
  case class FileHashResponse(hashResult: HashResult) extends BackendSpecificHasherResponse { override def hashes = Set(hashResult) }
  case class RuntimeAttributeHashesResponse(hashes: Set[HashResult]) extends BackendSpecificHasherResponse
}
