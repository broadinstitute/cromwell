package cromwell.services.metadata.impl.sns

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



class AwsSnsMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  implicit val ec: ExecutionContextExecutor = context.dispatcher

  //setup sns client
  val topicArn = "cromwell-metadata" //todo get this from serviceConfig

  //todo, should this be service config?
  val awsConfig: AwsConfiguration = AwsConfiguration(serviceConfig)
  val credentialsProviderChain: AwsCredentialsProviderChain =
    AwsCredentialsProviderChain.of(awsConfig.authsByName.values.map(_.provider()).toSeq :_*)

  lazy val snsClient: SnsClient = SnsClient.builder()
    .region(awsConfig.region.getOrElse(Region.US_EAST_1))
    .credentialsProvider(credentialsProviderChain)
    .build();

  def publishMessages(events: Iterable[MetadataEvent]): Future[Unit] = {
    import AwsSnsMetadataServiceActor.EnhancedMetadataEvents

    val eventsJson = events.toJson
    log.debug("Publishing to " + topicArn + ": " + eventsJson)

    def publish(): Unit = {
      snsClient.publish(PublishRequest.builder()
        .message(eventsJson.mkString("\n"))
        .topicArn(topicArn)
        .subject("cromwell-metadata-event")
        .build())
      () //return unit
    }

    Future(publish())
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


