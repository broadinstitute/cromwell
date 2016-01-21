package cromwell.logging

import akka.actor.{ActorContext, Actor, ActorLogging, Props}
import akka.util.Timeout
import org.slf4j.LoggerFactory
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
  val  eventLog = LoggerFactory.getLogger("cromwell.logging.BusinessLogEvent")

  override def info(message: String) = eventLog.info(message)
  override def info(message: String, t : Throwable) = eventLog.info(message,t)

  override def debug(message : String) = eventLog.debug(message)
  override def debug(message : String , t : Throwable) = eventLog.debug(message,t)

  override def warn(message: String) = eventLog.warn(message)
  override def warn(message: String, t : Throwable) = eventLog.warn(message,t)

  override def error(message: String) = eventLog.error(message)
  override def error(message: String , t : Throwable) = eventLog.error(message,t)

  override def trace(message: String) = eventLog.trace(message)
  override def trace(message: String, t : Throwable) = eventLog.trace(message,t)

  def logEvent(message: String): Unit = {
    info(message)
  }

  /**
    * This method receive events from different actor inside cromwell engine
    * like WorkflowActor , CallActor and so on.
    */
  val onBusinessLoggingEvent = (topic: WorkflowEvent, payload: Any) => Some(topic) collect {
    case event@(CallExecutionEvent | WorkflowExecutionEvent) =>
      logEvent(s"topic => $event and payload => $payload")
      info(s"<---------- received business event ------------>")
  }
}

sealed trait BusinessLoggingMessage
case object SubscribeToLogging extends BusinessLoggingMessage

/**
  * Specific logging implementation which is interested in listening
  * workflow event that could further categorize into workflow execution event
  * and call execution event.
  */
class BusinessLogging extends Actor with ActorLogging {
  this : BusinessLogEvent =>

  implicit val timeout = Timeout(10 seconds)
  implicit val actorSystem = context.system

  def actorRefFactory: ActorContext = context

  override def receive: Receive = {
    case SubscribeToLogging =>
      val subscriber = actorSystem.actorOf(Subscriber.props[WorkflowEvent](onBusinessLoggingEvent),"BusinessLogging")
      PubSubMediator.subscribe[WorkflowEvent](subscriber)
      log.info(s"received subscription message")
  }
}

object BusinessLogging {
  def props(): Props = Props(new BusinessLogging() with BusinessLogEvent)
}