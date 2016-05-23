package cromwell.services

import akka.actor.{ActorRef, Actor, ActorLogging}
import com.typesafe.config.Config
import cromwell.engine.db.DataAccess
import cromwell.services.MetadataServiceActor._
import wdl4s.values.{WdlArray, WdlValue}

import scala.concurrent.ExecutionContext.Implicits.global
import scala.language.postfixOps
import scala.util.{Failure, Success}

// TODO: PBE: Will not be MetadataServiceActor until circular dependencies fixed.
case class EngineMetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging {

  val dataAccess = DataAccess.globalDataAccess

  private def queryAndRespond(query: MetadataQuery) = {
    val sndr = sender()
    dataAccess.queryMetadataEvents(query) onComplete {
      case Success(m) => sndr ! MetadataLookupResponse(query, m)
      case Failure(t) => sndr ! MetadataServiceKeyLookupFailed(query, t)
    }
  }

  def receive = {
    case action@PutMetadataAction(event) =>
      val sndr = sender()
      dataAccess.addMetadataEvent(event) onComplete {
        case Success(_) => sndr ! MetadataPutAcknowledgement(action)
        case Failure(t) =>
          val msg = MetadataPutFailed(action, t)
          log.error(t, "Sending {} failure message {}", sndr, msg)
          sndr ! msg
      }
    case GetAllMetadataAction(workflowId) =>
      val query = MetadataQuery(workflowId, None, None)
      queryAndRespond(query)
    case GetMetadataQueryAction(query@MetadataQuery(_, _, _)) => queryAndRespond(query)
  }
}
