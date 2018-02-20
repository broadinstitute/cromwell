package cromwell.services.metadata.impl.splitpubsub

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata.impl.MetadataServiceActor
import cromwell.services.metadata.impl.WriteMetadataActor.CheckPendingWrites
import cromwell.services.metadata.impl.pubsub.PubSubMetadataServiceActor
import cromwell.util.GracefulShutdownHelper.ShutdownCommand


/**
  * A metadata service implementation which will route all put actions to pubsub, and all other actions to the standard MSA
  *
  * It is expected that the user will supply any desired config fields for both MetadataServiceActor and
  * PubSubMetadataServiceActor as the serviceConfig block will be passed along to both of them.
  */
class SplitPubSubMetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging {
  val standardMetadataActor: ActorRef = context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig))
  val pubSubMetadataActor: ActorRef = context.actorOf(PubSubMetadataServiceActor.props(serviceConfig, globalConfig))

  override def receive = {
    case action: PutMetadataAction =>
      pubSubMetadataActor forward action
    case action: PutMetadataActionAndRespond =>
      pubSubMetadataActor forward action

    case CheckPendingWrites =>
    case ShutdownCommand =>

    case anythingElse => standardMetadataActor forward anythingElse
  }
}
