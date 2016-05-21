package cromwell.services

import cromwell.core.WorkflowId
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

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

  /* TODO: PBE: No EngineMetadataServiceActor.props until circular dependencies fixed.
  def props(serviceConfig: Config, globalConfig: Config) = {
    Props(MetadataServiceActor(serviceConfig, globalConfig))
  }
  */
}
