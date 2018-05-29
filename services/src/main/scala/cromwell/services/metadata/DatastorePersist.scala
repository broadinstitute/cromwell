package cromwell.services.metadata


import akka.actor.{Actor, ActorLogging}
import cats.data.Kleisli
import cats.syntax.either._
import cats.effect.IO
import com.google.cloud.Timestamp
import com.google.cloud.datastore.{Datastore, DatastoreOptions, FullEntity, Key}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import net.ceedubs.ficus.Ficus._

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
  val pubSubConnection = datastore()

  override def receive = {
    case action: PutMetadataAction =>
      publishMessages(action.events).runAsync() unsafeRunAsync(_.fold({ e =>
        log.error(e, "Failed to post metadata: " + action.events)
      }, _ => ()))
    case action: PutMetadataActionAndRespond =>
      publishMessages(action.events).unsafeRunAsync({
        case Right(_) => action.replyTo ! MetadataWriteSuccess(action.events)
        case Left(e) => action.replyTo ! MetadataWriteFailure(e, action.events)
      })
  }

  protected[this] def datastore(): Datastore = DatastoreOptions.getDefaultInstance().getService();

  private def publishMessages(events: Iterable[MetadataEvent]): IO[Unit] = {

    val entities = events.map{
      case MetadataEvent(MetadataKey(workflowId, jobKey, key), Some(MetadataValue(value, valueType)), offsetDateTime) =>
        val taskKey = datastore.newKeyFactory.setKind("metadata").newKey(s"${workflowId.id.toString}${jobKey.map(jobKey => s"-$jobKey")}-$key")

        val instant = offsetDateTime.toInstant

        val timestamp = Timestamp.ofTimeSecondsAndNanos(instant.getEpochSecond,instant.getNano)

        FullEntity.
          newBuilder(taskKey).
          set("value", value).
          set("valueType", valueType.typeName).
          set("time", timestamp).build()
    }

    IO{pubSubConnection.put(entities.toSeq:_*)}.
      flatMap(_ => IO.shift(ec))
  }
}
