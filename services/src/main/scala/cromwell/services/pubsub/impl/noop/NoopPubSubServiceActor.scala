package cromwell.services.pubsub.impl.noop

import com.typesafe.config.Config
import cromwell.services.pubsub.PubSubServiceActor

/**
  * Intended for default cases which aren't actually publishing to any sort of external pubsub service, will just
  * drop any messages to the floor.
  */
class NoopPubSubServiceActor(val serviceConfig: Config, globalConfig: Config) extends PubSubServiceActor {
  override def publish(topic: String, message: String) = {
    logger.debug("Received request to publish " + message + " to topic " + topic + " but I choose not to run.")
  }
}
