package cromwell.services.metadata.impl.hybridpubsub

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata.impl.MetadataServiceActor
import cromwell.services.metadata.impl.pubsub.PubSubMetadataServiceActor

/**
  * A metadata service implementation which will function as a standard metadata service but also push all metadata
  * events to google pubsub.
  *
  * Under the hood it maintains its own MetadataServiceActor and PubSubMetadataServiceActor. All messages are routed
  * to the MetadataServiceActor. PutMetadataActions are also sent to the PubSubMetadataServiceActor. PutMetadataActionAndRespond
  * messages will be sent to the PubSubMetadataServiceActor as a standard PutMetadataAction, i.e. only the standard
  * metadata service will be ACKing the request.
  *
  * It is expected that the user will supply any desired config fields for both MetadataServiceActor and
  * PubSubMetadataServiceActor as the serviceConfig block will be passed along to both of them.
  */
class HybridPubSubMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends Actor
    with ActorLogging {
  val standardMetadataActor: ActorRef =
    context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))
  val pubSubMetadataActor: ActorRef =
    context.actorOf(PubSubMetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))

  override def receive = {
    case action: PutMetadataAction =>
      standardMetadataActor forward action
      pubSubMetadataActor forward action
    case action: PutMetadataActionAndRespond =>
      standardMetadataActor forward action
      pubSubMetadataActor forward PutMetadataAction(action.events)
    case anythingElse => standardMetadataActor forward anythingElse
  }
}
