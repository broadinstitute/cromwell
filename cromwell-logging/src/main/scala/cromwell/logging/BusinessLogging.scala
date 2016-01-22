package cromwell.logging

import akka.actor.{ActorContext, Actor, ActorLogging, Props}
import akka.util.Timeout
import org.slf4j.{Logger, LoggerFactory}
import scala.concurrent.duration._
import cromwell.pubsub.{Subscriber, PubSubMediator}
import scala.language.postfixOps

/**
  * Created by himanshu on 1/13/16.
  */
/**
  * LogWrapper is used to realize different types of logger
  * For example there is WorkflowLogger a specific implementation
  * of LogWrapper.
  */
trait LogWrapper {
  def debug(message : String)
  def debug(message : String , t : Throwable)
  def warn(message: String)
  def warn(message: String, t : Throwable)
  def error(message: String)
  def error(message: String , t : Throwable)
  def trace(message: String)
  def trace(message: String, t : Throwable)
  def info(message: String)
  def info(message: String, t : Throwable)
}

/**
  * Workflow execution event topics that business logging
  * subscribed to and expecting,
  */
sealed trait WorkflowEvent
case object CallExecutionEvent extends WorkflowEvent
case object WorkflowExecutionEvent extends WorkflowEvent

/**
  * To report Business logging received as event from
  * different components inside cromwell.
  */
trait BusinessLogEvent extends LogWrapper{
  def eventLogger() : Logger = LoggerFactory.getLogger("cromwell.logging.BusinessLogEvent")

  override def info(message: String) = eventLogger().info(message)
  override def info(message: String, t : Throwable) = eventLogger().info(message,t)

  override def debug(message : String) = eventLogger().debug(message)
  override def debug(message : String , t : Throwable) = eventLogger().debug(message,t)

  override def warn(message: String) = eventLogger().warn(message)
  override def warn(message: String, t : Throwable) = eventLogger().warn(message,t)

  override def error(message: String) = eventLogger().error(message)
  override def error(message: String , t : Throwable) = eventLogger().error(message,t)

  override def trace(message: String) = eventLogger().trace(message)
  override def trace(message: String, t : Throwable) = eventLogger().trace(message,t)

  def logEvent(message: String): Unit = {
    info(message)
  }

  /**
    * This method receive events from different actor inside cromwell engine
    * like WorkflowActor , CallActor and so on.
    */
  def onBusinessLoggingEvent() :(WorkflowEvent,Any) =>Option[Unit] = (topic: WorkflowEvent, payload: Any) => Some(topic) collect {
    case event@(CallExecutionEvent | WorkflowExecutionEvent) =>
      logEvent(s"topic => $event and payload => $payload")
  }
}

sealed trait BusinessLoggingMessage
case object SubscribeToLogging extends BusinessLoggingMessage
case object UnSubscribeToLogging extends BusinessLoggingMessage

/**
  * Specific logging implementation which is interested in listening
  * workflow event that could further categorize into workflow execution event
  * and call execution event.
  */
class BusinessLogging() extends Actor with ActorLogging {
  this : BusinessLogEvent with PubSubMediator=>

  implicit val timeout = Timeout(10 seconds)
  implicit val actorSystem = context.system

  def actorRefFactory: ActorContext = context
  val subscriber = actorSystem.actorOf(Subscriber.props[WorkflowEvent](onBusinessLoggingEvent()),s"BusinessLogging_${this.hashCode()}")

  override def receive: Receive = {
    case SubscribeToLogging =>
      subscribe[WorkflowEvent](subscriber)
      log.info(s"received subscription message")
    case UnSubscribeToLogging =>
      unsubscribe[WorkflowEvent](subscriber)
      log.info(s"received un subscription message")
  }
}

object BusinessLogging {
  def props(): Props = Props(new BusinessLogging() with BusinessLogEvent with PubSubMediator)
}