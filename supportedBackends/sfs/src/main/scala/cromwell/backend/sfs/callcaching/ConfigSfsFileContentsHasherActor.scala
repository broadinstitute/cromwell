package cromwell.backend.sfs.callcaching

import java.io.{File, FileInputStream}

import akka.actor.{Actor, ActorLogging, Props, Terminated}
import akka.routing.{ActorRefRoutee, RoundRobinRoutingLogic, Router}
import cromwell.backend.callcaching.BackendSpecificHasherActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.core.callcaching.{HashResult, HashValue, HashingFailedMessage}

import scala.util.{Failure, Success, Try}

class ConfigSfsFileContentsHasherActor extends Actor with ActorLogging {
  var router = {
    val routees = Vector.fill(5) {
      val r = context.actorOf(Props[SfsHashMakerWorker])
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
      val r = context.actorOf(Props[SfsHashMakerWorker])
      context watch r
      router = router.addRoutee(r)
  }
}

/**
  * Blocking worker. Warning: if this crashes its mailbox of hash requests will be lost!!
  */
private[callcaching] class SfsHashMakerWorker extends Actor with ActorLogging {
  override def receive = {
    case x: SingleFileHashRequest =>
       getMd5Result(x.file.valueString) match {
        case Success(md5Success) => sender ! FileHashResponse(HashResult(x.hashKey, HashValue(md5Success)))
        case Failure(t) => HashingFailedMessage(x.hashKey, t)
      }
  }

  private def getMd5Result(path: String) = {
    val md5Result = for {
      fileInputStream <- Try(new FileInputStream(new File(path)))
      // Use = here so that closeAndLog runs even if this fails:
      md5String = Try(org.apache.commons.codec.digest.DigestUtils.md5Hex(fileInputStream))
      _ = closeAndLog(fileInputStream)
    } yield md5String

    // Because we used = above, we must flatten the try down here
    md5Result.flatten
  }

  private def closeAndLog(fileInputStream: FileInputStream) = {
    Try(fileInputStream.close()) match {
      case Success(_) => // No need to log a success!
      case Failure(t) => log.error(s"Could not close file stream: $t")
    }
  }
}