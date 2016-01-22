package cromwell.pubsub

import akka.actor.{Props, Actor, ActorRef, ActorSystem}

/**
  * Created by himanshu on 1/12/16.
  */
/**
  * This is a wrapper to Eventstream as one of the event bus with SubChannel trait,
  * it is defined inside business logging component but later it needs to moved out
  * to separate place.
  */
//TODO:Needs to move some other place 
trait PubSubMediator{

  /**
    *Subscribes consumer function(topic,message) where topic is event structure and message is event data.
    * @param actorRef actor that needs to be subscribed,this Actor class is of type (T,Any) where
    *                 T is generic topic and Any is any payload
    * @param system system actor system for creating subscriber actor this is implicit
    * @tparam T
    * @return
    */
  def subscribe[T](actorRef: ActorRef)(implicit system:ActorSystem):Unit = {
    system.eventStream.subscribe(actorRef, classOf[(T, Any)])
  }

  /**
    * Producer produces event by publishing to Event Stream using publish method.
    * @param topic event topic this could be any t
    * @param payload
    * @param system
    * @tparam T
    */
  def publish[T](topic: T, payload: Any)(implicit system:ActorSystem):Unit = {
    system.eventStream.publish(topic, payload)
  }

  /**
    *This will un-subscribe from the akka system eventstream
    * @param system actor system implicit
    * @param actorRef actor that needs to be unsubscribe
    * @param topic optional topic that needs to be unsubscribe
    * @tparam T
    */
  def unsubscribe[T](actorRef: ActorRef,topic:Option[T] = None)(implicit system: ActorSystem): Unit = {
    topic match {
      case None => system.eventStream.unsubscribe(actorRef)
      case Some(x) => system.eventStream.unsubscribe(actorRef,classOf[(T,Any)])
    }
  }
}

object PubSubMediator extends PubSubMediator

/***
  * Subscriber actor allows consumer to subscribe for events, This actor is created by
  * parent actor how is intended to subscribe for events hence manage lifecycle of this actor.
  * @param f is a biFunction that receive workflow event as a topic and message as Any type.
  */
class Subscriber(f: (Any, Any) => Option[Unit]) extends Actor {
  override def receive = { case (topic: Any, payload: Any) => f(topic, payload) }
}

object Subscriber{
  def props[T](f :(T,Any) => Option[Unit]) : Props = Props(classOf[Subscriber],f)
}