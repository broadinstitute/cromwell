package cromwell.engine.workflow.lifecycle.execution.callcaching

import java.nio.file.Path
import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.JobKey
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.HashResultMessage
import cromwell.engine.workflow.lifecycle.execution.callcaching.FileHasherActor._

/**
  * How to bootstrap this TBD. Presumably a singleton.
  *
  * Contains a set of sub-actors, one per filesystem.
  *
  * Receives a set of files to hash, returns the hashes as a trickle.
  */
class FileHasherActor(filesystemSpecificActors: List[ActorRef]) extends Actor with ActorLogging {

  override def receive: Receive = {
    // Choose an appropriate hasher for each file and send it on
    case JobFileHashRequests(jobKey, files) =>
      files groupBy getAnAppropriateFilesystemSpecificActor foreach {
        case (actor, group) => actor.tell(JobFileHashRequests(jobKey, group), sender)
      }
    case cancellation: HashesNoLongerRequired =>
      filesystemSpecificActors foreach { _ ! cancellation }
  }

  // TODO: Implement me better (using the filesystemSpecificActors in the class constructor)
  // Should be QUICK (used in a receive block) otherwise re-engineer
  private val reallyStupidFileHasher = context.actorOf(ReallyStupidFileHasher.props)
  private def getAnAppropriateFilesystemSpecificActor(fileToHash: SingleFileHashRequest): ActorRef = reallyStupidFileHasher
}

object FileHasherActor {

  def props: Props = Props(new FileHasherActor(List.empty))

  sealed trait FileHasherCommand
  case class SingleFileHashRequest(hashKey: HashKey, file: Path)
  case class JobFileHashRequests(jobKey: JobKey, files: Iterable[SingleFileHashRequest]) extends FileHasherCommand
  case class HashesNoLongerRequired(jobKey: JobKey)

  case class FileHasherResponse(hashResult: HashResult) extends HashResultMessage { override def hashes = List(hashResult) }
}



@deprecated("This is really stupid", "always")
class ReallyStupidFileHasher extends Actor with ActorLogging {
  override def receive: Receive = {
    case JobFileHashRequests(jobKey, fileHashRequests) => {
      fileHashRequests foreach { hashRequest =>
        sender ! HashResult(hashRequest.hashKey, HashValue(hashRequest.file.toAbsolutePath.toString))
      }
    }
  }
}

@deprecated("This is really stupid", "always")
object ReallyStupidFileHasher {
  @deprecated("This is really stupid", "always")
  def props: Props = Props(new ReallyStupidFileHasher)
}
