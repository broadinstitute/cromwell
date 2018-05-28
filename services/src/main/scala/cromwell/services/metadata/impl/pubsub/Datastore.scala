package cromwell.services.metadata.impl.pubsub

import java.nio.file.{Files, Paths}

import akka.actor.{Actor, ActorLogging}
import cats.data.Validated.{Invalid, Valid}
import cats.instances.future._
import cats.syntax.functor._
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode.JsonFileFormat
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata._
import net.ceedubs.ficus.Ficus._
import org.broadinstitute.dsde.workbench.google.{GoogleCredentialModes, GooglePubSubDAO, HttpGooglePubSubDAO}
import com.google.cloud.datastore.{Datastore, DatastoreOptions, Entity}

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
class DatastoreMetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging {
  implicit val ec = context.dispatcher

  // The auth *must* be a service account auth but it might not be named service-account.
  val googleAuthName = serviceConfig.getOrElse("auth", "service-account")

  val googleProject = serviceConfig.as[String]("project")

  val pubSubAppName = serviceConfig.getOrElse("appName", "cromwell")
  val pubSubConnection = createPubSubConnection()

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

  protected[this] def createPubSubConnection(): Datastore = {
    /*
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
    */
        // Instantiates a client
    val datastore = DatastoreOptions.getDefaultInstance().getService();

    datastore
  }



  private def publishMessages(events: Iterable[MetadataEvent]): Future[Unit] = {
    import PubSubMetadataServiceActor.EnhancedMetadataEvents

    events.map{
      case MetadataEvent(MetadataKey(workflowId, jobKey, key), Some(MetadataValue(value, valueType)), offsetDateTime) =>
    }
    pubSubConnection.put()
  }
}



