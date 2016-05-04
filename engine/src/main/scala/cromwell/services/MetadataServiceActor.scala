package cromwell.services

import akka.actor.{Actor, Props}
import com.typesafe.config.Config
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.WorkflowId
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.services.MetadataServiceActor._
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import org.joda.time.DateTime

import scala.language.postfixOps

object MetadataServiceActor {

  final val MetadataServiceName = "MetadataService"

  case class MetadataJobKey(callFqn: String, index: ExecutionIndex, attempt: Int)
  object MetadataJobKey {
    def apply(backendJobKey: BackendJobDescriptorKey): MetadataJobKey = MetadataJobKey(backendJobKey.call.fullyQualifiedName, backendJobKey.index, backendJobKey.attempt)
  }

  case class MetadataKey(workflowId: WorkflowId, jobKey: Option[MetadataJobKey], key: String)
  case class MetadataValue(value: String)
  case class MetadataEvent(key: MetadataKey, value: MetadataValue, timestamp: DateTime = DateTime.now)

  case class MetadataQueryJobKey(callFqn: Option[String], index: Option[ExecutionIndex], attempt: Option[Int])
  object MetadataQueryJobKey {
    def forMetadataJobKey(jobKey: MetadataJobKey) = MetadataQueryJobKey(Option(jobKey.callFqn), Option(jobKey.index), Option(jobKey.attempt))
  }
  case class MetadataQuery(workflowId: Option[WorkflowId], jobKey: Option[MetadataQueryJobKey], key: Option[String])
  object MetadataQuery {
    def forWorkflow(workflowId: WorkflowId) = MetadataQuery(Option(workflowId), None, None)
    def forJob(workflowId: WorkflowId, jobKey: MetadataJobKey) = MetadataQuery(Option(workflowId), Option(MetadataQueryJobKey.forMetadataJobKey(jobKey)), None)
    def forKey(key: String) = MetadataQuery(None, None, Option(key))
  }

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
  def receive = {
    case PutMetadataAction(event) =>
      val result: MetadataPutAcknowledgement = ???
      sender ! result
    case GetAllMetadataAction(workflowId) =>
      val result: MetadataLookupResponse = ???
      sender ! result
    case GetMetadataQueryAction(MetadataQuery(workflowId, jobKey, key)) =>
      val result: MetadataLookupResponse = ???
      sender ! result
  }
}
