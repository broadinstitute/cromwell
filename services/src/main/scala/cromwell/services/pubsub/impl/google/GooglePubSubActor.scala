package cromwell.services.pubsub.impl.google

import java.security.PrivateKey
import java.util.{Base64, UUID}

import akka.NotUsed
import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.stream.{ActorMaterializer, FlowShape, OverflowStrategy}
import akka.stream.alpakka.googlecloud.pubsub.{PubSubMessage, PublishRequest}
import akka.stream.alpakka.googlecloud.pubsub.scaladsl.GooglePubSub
import akka.stream.scaladsl.{Flow, GraphDSL, Sink, Source, SourceQueueWithComplete, Unzip}
import better.files.File
import com.google.auth.oauth2.ServiceAccountCredentials
import cromwell.core.actor.StreamActorHelper
import cromwell.core.actor.StreamIntegration.StreamContext
import cromwell.services.pubsub.impl.google.GooglePubSubActor.{PubSubContext, PubSubRequest}

class FooPubSubActor extends GooglePubSubActor {
  val credentials = ServiceAccountCredentials.fromStream(new File("/Users/jgentry/broad-dsde-cromwell-dev-1ebbc164ce7e.json").newInputStream)

  override val privateKey = credentials.getPrivateKey
  override val clientEmail = credentials.getClientEmail
  override val projectId = credentials.getProjectId
  override val topic = "test-ps"
  override val apiKey = "some bullshit"

  override val sourceQueueSize = 100
}

/**
  * FIXME: Move to cloudSupport
  *
  * Provides capability to receive messages to publish to a Google Cloud PubSub topic.
  *
  * FIXME: What if topic has not been created?
  */
trait GooglePubSubActor extends Actor with ActorLogging with StreamActorHelper[PubSubContext] {
  implicit val actorSystem = context.system
  implicit val materializer = ActorMaterializer()
  implicit val ec = actorSystem.dispatcher

  def sourceQueueSize: Int

  // Fields tied to the Alpakka connector
  def projectId: String
  def apiKey: String
  def privateKey: PrivateKey
  def clientEmail: String
  def topic: String
  def parallelism: Int = 1

  override protected def actorReceive: Receive = {
    case request: PubSubRequest => println("I RECEIVED: " + request)//sendToStream(PubSubContext(request, sender()))
  }

  protected val publishSource: Source[PubSubContext, SourceQueueWithComplete[PubSubContext]] = {
    Source.queue[PubSubContext](sourceQueueSize, OverflowStrategy.dropNew)
  }

  val publishFlow: Flow[PublishRequest, Seq[String], NotUsed] = GooglePubSub.publish(projectId, apiKey, clientEmail, privateKey, topic)

  val someFlow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    val input = builder.add(Flow[(PubSubRequest, PubSubContext)])

    val unzipper = builder.add(Unzip[PubSubRequest, PubSubContext])
    val toPublishRequest = unzipper.out0 map { p => PublishRequest(List(PubSubMessage(new String(Base64.getEncoder.encode(p.message.getBytes)), UUID.randomUUID().toString))) }

    val requestsFlow = builder.add(publishFlow)

    input ~> unzipper.in
             toPublishRequest.outlet ~> requestsFlow ~> Sink.ignore

    FlowShape[(PubSubRequest, PubSubContext), PubSubContext](input.in, unzipper.out1)
  }

  override protected def streamSource = publishSource map { context => (context.request, context) }

  override val replySink = Sink.foreach[(Any, PubSubContext)] {
    case (_, _) => () // response and commandContext
  }
}

object GooglePubSubActor {
  final case class PubSubRequest(message: String) extends AnyVal

  case class PubSubContext(request: PubSubRequest, replyTo: ActorRef) extends StreamContext
}