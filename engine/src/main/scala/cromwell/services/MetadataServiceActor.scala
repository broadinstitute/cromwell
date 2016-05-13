package cromwell.services

import akka.actor.{Actor, Props}
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.engine.db.DataAccess
import cromwell.services.MetadataServiceActor._
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

import scala.concurrent.ExecutionContext.Implicits.global
import scala.language.postfixOps
import scala.util.{Failure, Success}


object MetadataServiceActor {

  final val MetadataServiceName = "MetadataService"

  trait MetadataServiceMessage
  /**
    * Command Actions
    */
  trait MetadataServiceAction extends MetadataServiceMessage with ServiceRegistryMessage {
    def serviceName = MetadataServiceName
  }
  case class PutMetadataAction(event: MetadataEvent) extends MetadataServiceAction
  case class GetAllMetadataAction(workflowId: WorkflowId) extends MetadataServiceAction
  case class GetMetadataQueryAction(key: MetadataQuery) extends MetadataServiceAction

  /**
    * Responses
    */
  trait MetadataServiceResponse extends MetadataServiceMessage
  case class MetadataPutAcknowledgement(putRequest: PutMetadataAction) extends MetadataServiceResponse
  case class MetadataPutFailed(putRequest: PutMetadataAction, reason: Throwable) extends MetadataServiceResponse

  case class MetadataLookupResponse(query: MetadataQuery, eventList: Seq[MetadataEvent]) extends MetadataServiceResponse
  case class MetadataServiceKeyLookupFailed(query: MetadataQuery, reason: Throwable) extends MetadataServiceResponse

  def props(serviceConfig: Config, globalConfig: Config) = {
    Props(MetadataServiceActor(serviceConfig, globalConfig))
  }
}

case class MetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor {

  private val dataAccess = DataAccess.globalDataAccess

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
        case Failure(t) => sndr ! MetadataPutFailed(action, t)
      }
    case GetAllMetadataAction(workflowId) =>
      val query = MetadataQuery(workflowId, None, None)
      queryAndRespond(query)
    case GetMetadataQueryAction(query@MetadataQuery(_, _, _)) => queryAndRespond(query)
  }
}
