package cromwell.services

import akka.actor.ActorRef
import cromwell.core.WorkflowId
import cromwell.services.MetadataServiceActor.PutMetadataAction
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import wdl4s.values._

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

object MetadataServiceActorImplicits {
  implicit class EnhancedServiceRegistryActorForMetadata(val actor: ActorRef) extends AnyVal {
    def pushWdlValueMetadata(metadataKey: MetadataKey, output: WdlValue): Unit = output match {
      case WdlArray(_, valueSeq) =>
        val zippedSeq = valueSeq.zipWithIndex
        zippedSeq foreach { case (value, index) => actor.pushWdlValueMetadata(metadataKey.copy(key = s"${metadataKey.key}[$index]"), value) }
      case WdlMap(_, valueMap) =>
        valueMap foreach { case (key, value) => actor.pushWdlValueMetadata(metadataKey.copy(key = metadataKey.key + s":${key.valueString}"), value) }
      case value =>
        actor ! PutMetadataAction(MetadataEvent(metadataKey, MetadataValue(value)))
    }

    def pushThrowableMetadata(metadataKey: MetadataKey, t: Throwable): Unit = {
      actor ! PutMetadataAction(MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:message"), MetadataValue(t.getMessage)))
      Option(t.getCause) foreach { cause => pushThrowableMetadata(metadataKey.copy(key = s"${metadataKey.key}:causedBy"), cause) }
    }
  }
}
