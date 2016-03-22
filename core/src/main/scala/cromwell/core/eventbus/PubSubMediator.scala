package cromwell.core.eventbus

import akka.actor.{ActorRef, ActorSystem}

final case class EventMsg[A](topic: A, payload: Any)

/**
  * This is a wrapper to EventStream as one of the event bus with SubChannel trait.
  */
trait PubSubMediator {
  /**
    * Subscribes consumer function(topic, message) where topic is event structure and message is event data.
    *
    * @param actorRef actor that needs to be subscribed,this Actor class is of type (A, Any) where
    *                 A is generic topic and Any is any payload
    * @param system   system actor system for creating subscriber actor this is implicit
    * @tparam A
    * @return
    */
  def subscribe[A](actorRef: ActorRef)(implicit system: ActorSystem): Unit = {
    system.eventStream.subscribe(actorRef, classOf[EventMsg[A]])
  }

  /**
    * Producer produces event by publishing to Event Stream using publish method.
    *
    * @param system
    * @tparam A
    */
  def publish[A](eventMsg: EventMsg[A])(implicit system: ActorSystem): Unit = {
    system.eventStream.publish(eventMsg)
  }

  /**
    * This will un-subscribe from the akka system EventStream
    *
    * @param system   actor system implicit
    * @param actorRef actor that needs to be unsubscribe
    * @param topic    optional topic that needs to be unsubscribe
    * @tparam A
    */
  def unsubscribe[A](actorRef: ActorRef, topic: Option[A] = None)(implicit system: ActorSystem): Unit = {
    topic match {
      case Some(x) => system.eventStream.unsubscribe(actorRef, classOf[EventMsg[A]])
      case None => system.eventStream.unsubscribe(actorRef)
    }
  }
}