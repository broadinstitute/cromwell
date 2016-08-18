package cromwell.backend.impl.jes.callcaching

import akka.actor.{Actor, ActorLogging, Props, Terminated}
import akka.routing.{ActorRefRoutee, RoundRobinRoutingLogic, Router}
import cromwell.backend.callcaching.BackendSpecificHasherActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.backend.impl.jes.JesBackendInitializationData
import cromwell.core.callcaching.{HashResult, HashValue, HashingFailedMessage}

import scala.util.{Failure, Success, Try}

class JesFileContentsHasherActor extends Actor with ActorLogging {
  var router = {
    val routees = Vector.fill(100) {
      val r = context.actorOf(Props[GcsHashFetcherWorker])
      context watch r
      ActorRefRoutee(r)
    }
    Router(RoundRobinRoutingLogic(), routees)
  }

  def receive = {
    case w: SingleFileHashRequest =>
      router.route(w, sender())
    case Terminated(a) =>
      router = router.removeRoutee(a)
      val r = context.actorOf(Props[GcsHashFetcherWorker])
      context watch r
      router = router.addRoutee(r)
  }
}

/**
  * Blocking worker. Warning! If this actor dies then its mailbox of hash requests will be lost
  */
private[callcaching] class GcsHashFetcherWorker extends Actor with ActorLogging {
  override def receive = {
    case x: SingleFileHashRequest =>

      // Create the path with the filesystem in the initialization data:
      getCrc32c(x) match {
        case Success(crc32cSuccess) => sender ! FileHashResponse(HashResult(x.hashKey, HashValue(crc32cSuccess)))
        case Failure(t) => HashingFailedMessage(x.hashKey, t)
      }
  }

  private def getCrc32c(singleFileHashRequest: SingleFileHashRequest): Try[String] = {
    def usingJesInitData(jesInitData: JesBackendInitializationData) = for {
      path <- Try(jesInitData.workflowPaths.fileSystemWithGenomicsAuth.getPath(singleFileHashRequest.file.valueString))
      crc32c <- Try(jesInitData.workflowPaths.fileSystemProvider.crc32cHash(path))
    } yield crc32c

    singleFileHashRequest.initializationData match {
      case Some(jesInitData: JesBackendInitializationData) => usingJesInitData(jesInitData)
      case _ => Failure(new IllegalArgumentException("Need JesBackendInitializationData to generate a GCS CRC32C hash"))
    }
  }
}