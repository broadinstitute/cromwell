package cromwell.services.metadata.impl.aws

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata.impl.MetadataServiceActor


/**
  * A metadata service implementation which will function as a standard metadata service but also push all metadata
  * events to AWS EventBridge. This class closely follows the pattern established in the
  * HybridPubSubMetadataServiceActor
  *
  * Under the hood it maintains its own MetadataServiceActor and AwsEventBridgeMetadataServiceActor. All messages are routed
  * to the MetadataServiceActor. PutMetadataActions are also sent to the AwsEventBridgeMetadataServiceActor. PutMetadataActionAndRespond
  * messages will be sent to the EventBridgeMetadataServiceActor as a standard PutMetadataAction, i.e. only the standard
  * metadata service will be ACKing the request.
  *
  * To use this actor something similar to the following should be present in the cromwell.conf file:
  * <pre>
  * services {
  *   MetadataService {
  *     class="cromwell.services.metadata.impl.aws.HybridEventBridgeMetadataServiceActor"
  *     config {
  *       aws {
  *         application-name = "cromwell"
  *         auths = [{
  *           name = "default"
  *           scheme = "default"
  *         }]
  *         region = "us-east-1"
  *         busName = "cromwell-metadata"
  *       }
  *     }
  *   }
  * }
  * </pre>
  *
  * @see cromwell.services.metadata.impl.aws.AwsEventBridgeMetadataServiceActor
  */
class HybridEventBridgeMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging  {
  val standardMetadataActor: ActorRef = context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))
  val awsEventBridgeMetadataActor: ActorRef = context.actorOf(AwsEventBridgeMetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))

  override def receive = {
    case action: PutMetadataAction =>
      standardMetadataActor forward action
      awsEventBridgeMetadataActor forward action
    case action: PutMetadataActionAndRespond =>
      standardMetadataActor forward action
      awsEventBridgeMetadataActor forward PutMetadataAction(action.events)
    case anythingElse => standardMetadataActor forward anythingElse
  }
}

