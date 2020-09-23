package cromwell.services.metadata.impl.pubsub

import java.nio.file.{Files, Paths}

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.functor._
import cats.instances.future._
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode.JsonFileFormat
import cromwell.core.Dispatcher._
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata._
import net.ceedubs.ficus.Ficus._
import org.broadinstitute.dsde.workbench.google.{GoogleCredentialModes, GooglePubSubDAO, HttpGooglePubSubDAO}
import spray.json._

import scala.concurrent.Future
import scala.util.{Failure, Success}

/**
  * A *write-only* metadata service implementation which pushes all events to Google PubSub. The expectation is that
  * metadata reads are being handled outside of this Cromwell instance.
  *
  * For now it requires a google service account auth with a PEM file. If the underlying PubSub library is updated to
  * allow for other forms of auth this becomes viable for us. Note that currently Cromwell complains that PEM is
  * deprecated, so we might need to contribute a patch to workbench-lib's pubsub tooling.
  *
  * This is intentionally not documented, considering there exists no way to consume/read metadata if this is active
  * anyone using this is by definition more advanced than the Cromwell developers
  *
  * Potential TODOs to consider if moving this forward from proof-of-concept:
  *  - This will probably need some measure of batching, need to see it in action before proceeding down that path
  *  - Graceful shutdown could also be useful. The pubsub calls will autoretry in the face of transient errors which
  *     means there might be stuff hanging around on a shutdown request.
  *  - Currently service actors which fail to instantiate properly don't bring down Cromwell. It's easier to screw
  *     it up here, we might want to revisit that strategy
  *  - Revisit the Alpakka connector instead of the workbench-lib version. There was a bug they fixed but hadn't released
  *     it yet. Also I'm pretty sure that same bug was repeated in a few other places.
  */
class PubSubMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  implicit val ec = context.dispatcher

  // The auth *must* be a service account auth but it might not be named service-account.
  val googleAuthName = serviceConfig.getOrElse("auth", "service-account")

  val googleProject = serviceConfig.as[String]("project")
  val pubSubTopicName = serviceConfig.getOrElse("topic", "cromwell-metadata")
  val pubSubSubscriptionName = serviceConfig.getAs[String]("subscription")
  val pubSubAppName = serviceConfig.getOrElse("appName", "cromwell")

  val pubSubConnection = createPubSubConnection()
  createTopicAndSubscription().failed foreach { e => throw e }

  override def receive = {
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

  protected[this] def createPubSubConnection(): GooglePubSubDAO = {
    implicit val as = context.system

    val googleConfig = GoogleConfiguration(globalConfig)

    // This class requires a service account auth due to the library used
    val googleAuth = googleConfig.auth(googleAuthName) match {
      case Valid(a: ServiceAccountMode) => a
      case Valid(doh) => throw new IllegalArgumentException(s"Unable to configure PubSubMetadataServiceActor: ${doh.name} was not a service account auth")
      case Invalid(e) => throw new IllegalArgumentException("Unable to configure PubSubMetadataServiceActor: " + e.toList.mkString(", "))
    }
    val jsonAuth = googleAuth.fileFormat match {
      case j: JsonFileFormat =>
        val jsonString = new String(Files.readAllBytes(Paths.get(j.file)))
        GoogleCredentialModes.Json(jsonString)
      case _ => throw new IllegalArgumentException("Unable to configure PubSubMetadataServiceActor: the service account must supply a JSON file")
    }

    // The metric name is pretty useless. It's an artifact of the workbench libraries which we're not using
    new HttpGooglePubSubDAO(pubSubAppName, jsonAuth, "cromwell.pubsub", googleProject)
  }

  /**
    * Creates the specified topic, and optionally the subscription if it was defined. We are ignoring the bool return
    * type from the HttpGooglePubSubDAO return types as it is not considered an error if the topic/subscription already
    * exist (in fact, this is a likely circumstance).
    */
  private def createTopicAndSubscription(): Future[Unit] = {
    log.info("Ensuring topic " + pubSubTopicName + " exists")
    pubSubConnection.createTopic(pubSubTopicName).flatMap(_ => createSubscription()).void
  }

  /**
    * Creates a subscription if pubsubTopicName is defined. Note that we're ignoring the bool in the Future[bool]
    * return type of HttpGooglePubSubDAO.createSubscription and only caring about the success of the Future. It is not
    * considered an error if the subscription already exists (in fact, this is a likely circumstance).
    */
  private def createSubscription(): Future[Unit] = {
    pubSubSubscriptionName match {
      case Some(name) =>
        log.info("Creating subscription " + name)
        pubSubConnection.createSubscription(pubSubTopicName, name).void
      case None =>
        log.info("Not creating a subscription")
        Future.successful(())
    }
  }

  private def publishMessages(events: Iterable[MetadataEvent]): Future[Unit] = {
    import PubSubMetadataServiceActor.EnhancedMetadataEvents

    val eventsJson = events.toJson
    log.debug("Publishing to " + pubSubTopicName + ": " + eventsJson)
    pubSubConnection.publishMessages(pubSubTopicName, eventsJson)
  }
}

object PubSubMetadataServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = {
    Props(new PubSubMetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
  }

  implicit class EnhancedMetadataEvents(val e: Iterable[MetadataEvent]) extends AnyVal {
    import MetadataJsonSupport._

    def toJson: Seq[String] = e.map(_.toJson.toString()).toSeq
  }
}

