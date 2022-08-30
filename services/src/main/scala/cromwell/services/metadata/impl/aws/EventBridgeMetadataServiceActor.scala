package cromwell.services.metadata.impl.aws

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.cloudsupport.aws.AwsConfiguration
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import software.amazon.awssdk.auth.credentials.AwsCredentialsProviderChain
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.eventbridge.EventBridgeClient
import software.amazon.awssdk.services.eventbridge.model.PutEventsRequest
import software.amazon.awssdk.services.eventbridge.model.PutEventsRequestEntry
import spray.json.enrichAny

import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}


/**
  * An actor that publishes metadata events to AWS EventBridge
  * @param serviceConfig the source of service config information
  * @param globalConfig the source of global config information
  * @param serviceRegistryActor the actor for registering services
  * @see cromwell.services.metadata.impl.aws.HybridEventBridgeMetadataServiceActor
  */
class AwsEventBridgeMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  implicit val ec: ExecutionContextExecutor = context.dispatcher

  //setup EB client
  val busName: String = serviceConfig.getString("aws.busName")

  val awsConfig: AwsConfiguration = AwsConfiguration(globalConfig)
  val credentialsProviderChain: AwsCredentialsProviderChain =
    AwsCredentialsProviderChain.of(awsConfig.authsByName.values.map(_.provider()).toSeq :_*)

  lazy val eventBrClient : EventBridgeClient = EventBridgeClient.builder()
    .region(awsConfig.region.getOrElse(Region.US_EAST_1))
    .credentialsProvider(credentialsProviderChain)
    .build();

  def publishMessages(events: Iterable[MetadataEvent]): Future[Unit] = {
    import AwsEventBridgeMetadataServiceActor.EnhancedMetadataEvents

    val eventsJson = events.toJson
    //if there are no events then don't publish anything
    if( eventsJson.length < 1) { return Future(())}
    log.debug(f"Publishing to $busName : $eventsJson")

    val reqEntry = PutEventsRequestEntry.builder()
      .eventBusName(busName)
      .source("cromwell")
      .detailType("cromwell-metadata-event")
      .detail(eventsJson.mkString(","))
      .build()

    val eventsRequest = PutEventsRequest.builder()
      .entries(reqEntry)
      .build()

    Future {
      eventBrClient.putEvents(eventsRequest)
      () //return unit
    }
  }

  override def receive: Receive = {
    case action: PutMetadataAction =>
      publishMessages(action.events).failed foreach { e =>
        log.error(e, "Failed to post metadata: " + action.events)
      }
    case action: PutMetadataActionAndRespond =>
      publishMessages(action.events) onComplete {
        case Success(_) => action.replyTo ! MetadataWriteSuccess(action.events)
        case Failure(e) => action.replyTo ! MetadataWriteFailure(e, action.events)
      }
  }
}

object AwsEventBridgeMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props = {
    Props(new AwsEventBridgeMetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
  }

  implicit class EnhancedMetadataEvents(val e: Iterable[MetadataEvent]) extends AnyVal {
    import cromwell.services.metadata.MetadataJsonSupport._

    def toJson: Seq[String] = e.map(_.toJson.toString()).toSeq
  }
}

