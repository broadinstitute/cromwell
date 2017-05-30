package cromwell.core.actor

import akka.actor.{Actor, ActorLogging}
import akka.pattern.pipe
import akka.stream.{ActorAttributes, ActorMaterializer, QueueOfferResult, Supervision}
import akka.stream.QueueOfferResult.{Dropped, Enqueued, QueueClosed}
import akka.stream.scaladsl.{Sink, Source, SourceQueueWithComplete}
import cromwell.core.actor.StreamActorHelper.{ActorRestartException, StreamCompleted, StreamFailed}
import cromwell.core.actor.StreamIntegration._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object StreamActorHelper {
  private [actor] case class StreamFailed(failure: Throwable)
  private [actor] case object StreamCompleted
  class ActorRestartException(throwable: Throwable) extends RuntimeException(throwable)
}

trait StreamActorHelper[T <: StreamContext] { this: Actor with ActorLogging =>

  implicit def ec: ExecutionContext

  implicit def materializer: ActorMaterializer

  private val decider: Supervision.Decider = _ => Supervision.Resume
  
  private val replySink = Sink.foreach[(Any, T)] {
    case (response, commandContext) =>
      val reply = commandContext.clientContext map { (_, response) } getOrElse response
      commandContext.replyTo ! reply
  }

  protected def actorReceive: Receive
  
  protected def streamSource: Source[(Any, T), SourceQueueWithComplete[T]]

  override def receive = streamReceive.orElse(actorReceive)

  private [actor] lazy val stream = {
    streamSource
      .to(replySink)
      .withAttributes(ActorAttributes.supervisionStrategy(decider))
      .run()
  }

  override def preStart(): Unit = {
    stream.watchCompletion() onComplete {
      case Success(_) =>
        self ! StreamCompleted
      case Failure(failure) =>
        self ! StreamFailed(failure)
    }
  }

  def sendToStream(commandContext: T) = {
    val enqueue = stream offer commandContext map {
      case Enqueued => EnqueueResponse(Enqueued, commandContext)
      case other => EnqueueResponse(other, commandContext)
    } recoverWith {
      case t => Future.successful(FailedToEnqueue(t, commandContext))
    }

    pipe(enqueue) to self
    ()
  }
  
  private def backpressure(commandContext: StreamContext) = {
    val originalRequest = commandContext.clientContext map { _ -> commandContext.request } getOrElse commandContext.request
    commandContext.replyTo ! BackPressure(originalRequest)
  }

  private def streamReceive: Receive = {
    case EnqueueResponse(Enqueued, commandContext: T @unchecked) => // Good !
    case EnqueueResponse(Dropped, commandContext) => backpressure(commandContext)
      
      // In any of the cases below, the stream is in a failed state, which will be caught by the watchCompletion hook and the
      // actor will be restarted
    case EnqueueResponse(QueueClosed, commandContext) => backpressure(commandContext)
    case EnqueueResponse(QueueOfferResult.Failure(failure), commandContext) => backpressure(commandContext)
    case FailedToEnqueue(throwable, commandContext) => backpressure(commandContext)
      
      // Those 2 cases should never happen, as long as the strategy is Resume, but in case it does...
    case StreamCompleted => restart(new IllegalStateException("Stream was completed unexpectedly"))
    case StreamFailed(failure) => restart(failure)
  }

  /** Throw the exception to force the actor to restart so it can be back in business
    * IMPORTANT: Make sure the supervision strategy for this actor is Restart
    */
  private def restart(throwable: Throwable) = {
    throw new ActorRestartException(throwable)
  }
}
