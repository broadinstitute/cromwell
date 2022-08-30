package cromwell.services.metadata.impl.aws

import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.cloudsupport.aws.AwsConfiguration
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import software.amazon.awssdk.auth.credentials.AwsCredentialsProviderChain
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.sns.SnsClient
import software.amazon.awssdk.services.sns.model.PublishRequest
import spray.json.enrichAny

import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}


/**
  * An actor that publishes metadata events to AWS SNS
  * @param serviceConfig the source of service config information
  * @param globalConfig the source of global config information
  * @param serviceRegistryActor the actor for registering services
  * @see cromwell.services.metadata.impl.sns.HybridSnsMetadataServiceActor
  */
class AwsSnsMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  implicit val ec: ExecutionContextExecutor = context.dispatcher

  //setup sns client
  val topicArn: String = serviceConfig.getString("aws.topicArn")

  val awsConfig: AwsConfiguration = AwsConfiguration(globalConfig)
  val credentialsProviderChain: AwsCredentialsProviderChain =
    AwsCredentialsProviderChain.of(awsConfig.authsByName.values.map(_.provider()).toSeq :_*)

  lazy val snsClient: SnsClient = SnsClient.builder()
    .region(awsConfig.region.getOrElse(Region.US_EAST_1))
    .credentialsProvider(credentialsProviderChain)
    .build()

  def publishMessages(events: Iterable[MetadataEvent]): Future[Unit] = {
    import AwsSnsMetadataServiceActor.EnhancedMetadataEvents

    val eventsJson = events.toJson
    //if there are no events then don't publish anything
    if( eventsJson.length < 1) { return Future(())}
    log.debug(f"Publishing to $topicArn : $eventsJson")

    val message = PublishRequest.builder()
      .message("[" + eventsJson.mkString(",") + "]")
      .topicArn(topicArn)
      .subject("cromwell-metadata-event")

    if (topicArn.endsWith(".fifo")) {
      message
        .messageGroupId("cromwell")
        .messageDeduplicationId(UUID.randomUUID().toString())
    }

    Future {
      snsClient.publish(message
        .build())
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

object AwsSnsMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef): Props = {
    Props(new AwsSnsMetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
  }

  implicit class EnhancedMetadataEvents(val e: Iterable[MetadataEvent]) extends AnyVal {
    import cromwell.services.metadata.MetadataJsonSupport._

    def toJson: Seq[String] = e.map(_.toJson.toString()).toSeq
  }
}

