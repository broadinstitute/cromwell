package cromwell.services.pubsub.impl.google

import akka.stream.ActorMaterializer
import akka.stream.alpakka.googlecloud.pubsub.PubSubMessage
import akka.stream.alpakka.googlecloud.pubsub.scaladsl.GooglePubSub
import com.typesafe.config.Config
import cromwell.services.pubsub.PubSubServiceActor

class GooglePubSubServiceActor(val serviceConfig: Config, globalConfig: Config) extends PubSubServiceActor {
  implicit val actorSystem = context.system
  implicit val mat = ActorMaterializer()

  val z = PubSubMessage(messageId = ???, data = ???)
  val x = GooglePubSub.publish
}
