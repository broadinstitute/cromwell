package cromwell.pubsub


/**
  * Created by himanshu on 1/12/16.
  */
import akka.actor.{Actor, ActorSystem, Props}

import scala.reflect.ClassTag

/***
  * Subscriber actor allows consumer to subscribe for events, This actor is created by
  * parent actor how is intended to subscribe for events hence manage lifecycle of this actor.
  * @param f is a biFunction that receive workflow event as a topic and message as Any type.
  */
sealed class Subscriber[T : ClassTag](f: (T, Any) => Unit) extends Actor {
  override def receive = { case (topic: T, payload: Any) => f(topic, payload) }
}

object EventStream{

  /**
    * Subscribes consumer function(topic,message) where topic is event structure and message is event data.
    * @param system actor system for creating subscriber actor
    * @param f it is a function of parameters topic of type T that can be defined and message of type Any
    * @param name is subscriber name
    * @return
    */
  def subscribe[T : ClassTag](system:ActorSystem,f: (T, Any) => Option[Unit], name: String) = {
    val props = Props(classOf[Subscriber[T]], f)
    val subscriber = system.actorOf(props, name = name)
    system.eventStream.subscribe(subscriber, classOf[(T, Any)])
  }

  /**
    * Producer produces event by publishing to Event Stream using publish method.
    * @param system actor system for publishing to Event Stream
    * @param topic event name
    * @param payload event data
    * @tparam T is event type
    */
  def publish[T](system:ActorSystem,topic: T, payload: Any) {
    system.eventStream.publish(topic, payload)
  }

}