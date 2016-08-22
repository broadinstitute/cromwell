package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor._
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.DockerHashLookupWorkerActor.{DockerHashLookupCommand, DockerHashLookupKey, DockerHashLookupResponse}

import scala.util.{Failure, Success, Try}

/**
  * Blocking receive method. Recommended use is within a router.
  */
class DockerHashLookupWorkerActor extends Actor with ActorLogging {
  override def receive = {
    case x: DockerHashLookupCommand =>
      getDockerHash(x) match {
        case Success(dockerHashLookupSuccess) => sender ! DockerHashLookupResponse(HashResult(DockerHashLookupKey, HashValue(dockerHashLookupSuccess)))
        case Failure(t) => sender ! HashingFailedMessage(DockerHashLookupKey, t)
      }
  }

  def getDockerHash(lookupCommand: DockerHashLookupCommand): Try[String] = Failure(new NotImplementedError("No docker hash lookup yet"))
}

object DockerHashLookupWorkerActor {
  object DockerHashLookupKey extends HashKey("runtime attribute: docker (HASH)")
  case class DockerHashLookupCommand(name: String)
  case class DockerHashLookupResponse(hashResult: HashResult) extends SuccessfulHashResultMessage { override val hashes = Set(hashResult) }
}
