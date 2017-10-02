package cromwell.services.pubsub.impl.google

import java.security.PrivateKey

import akka.NotUsed
import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.stream.alpakka.googlecloud.pubsub
import akka.stream.{ActorAttributes, ActorMaterializer, OverflowStrategy}
import akka.stream.alpakka.googlecloud.pubsub.PublishRequest
import akka.stream.alpakka.googlecloud.pubsub.scaladsl.GooglePubSub
import akka.stream.scaladsl.SourceQueueWithComplete
import akka.stream.scaladsl.{Flow, Source}
import cromwell.core.Dispatcher
import cromwell.core.actor.StreamActorHelper
import cromwell.core.actor.StreamIntegration.StreamContext
import cromwell.services.pubsub.impl.google.GooglePubSubActor.{PubSubContext, PubSubRequest}

/**
  * Provides capability to receive messages to publish to a Google Cloud PubSub topic.
  *
  * FIXME: What if topic has not been created?
  */
trait GooglePubSubActor extends Actor with ActorLogging with StreamActorHelper[PubSubContext] {
  implicit val actorSystem = context.system
  implicit val mat = ActorMaterializer()

  def sourceQueueSize: Int

  // Fields tied to the Alpakka connector
  def projectId: String
  def apiKey: String
  def privateKey: PrivateKey
  def clientEmail: String
  def topic: String
  def parallelism: Int = 1

  override protected def actorReceive: Receive = {
    case request: PubSubRequest => sendToStream(PubSubContext(request, sender()))
  }

  protected val publishSource: Source[PubSubRequest, SourceQueueWithComplete[PubSubRequest]] = {
    Source.queue[PubSubRequest](sourceQueueSize, OverflowStrategy.dropNew)
  }

  val publishFlow: Flow[PublishRequest, Seq[String], NotUsed] = GooglePubSub.publish(projectId, apiKey, clientEmail, privateKey, topic)

  val z = publishSource.map(x => PublishRequest(List(pubsub.PubSubMessage(???, ???)))).via(publishFlow)

  override protected def streamSource = z
}

object GooglePubSubActor {
  final case class PubSubRequest(message: String) extends AnyVal

  case class PubSubContext(request: PubSubRequest, replyTo: ActorRef) extends StreamContext
}