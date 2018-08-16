package cromwell.engine.workflow.tokens

import akka.actor.ActorRef
import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.TokenQueue.{DequeueResult, LeasedActor, ShuffledQueueAndLease, TokenQueuePlaceholder}
import cromwell.engine.workflow.tokens.UnhoggableTokenPool.UnhoggableTokenPoolResult
import io.github.andrebeat.pool.Lease

import scala.collection.immutable.Queue

object TokenQueue {
  case class DequeueResult(leasedActor: Option[LeasedActor], tokenQueue: TokenQueue)
  case class LeasedActor(actor: ActorRef, lease: Lease[JobExecutionToken])
  def apply(tokenType: JobExecutionTokenType) = new TokenQueue(Queue.empty, UnhoggableTokenPool(tokenType))
  final case class TokenQueuePlaceholder(actor: ActorRef, hogGroup: String)
  final case class ShuffledQueueAndLease(lease: Option[Lease[JobExecutionToken]], newQueue: Queue[TokenQueuePlaceholder])
}

/**
  * A queue assigned to a pool.
  * Elements can be dequeued if the queue is not empty and there are tokens in the pool
  */
final case class TokenQueue(queue: Queue[TokenQueuePlaceholder], private [tokens] val pool: UnhoggableTokenPool) {
  val tokenType = pool.tokenType

  /**
    * Size of the queue
    */
  def size = queue.size

  /**
    * Enqueues an actor
    *
    * @return the new token queue
    */
  def enqueue(placeholder: TokenQueuePlaceholder): TokenQueue = {
    copy(queue = queue.enqueue(placeholder))
  }

  /**
    * Returns a dequeue'd actor if one exists and there's a token available for it
    * Returns an updated token queue based on this request (successful queuees get removed, hogs get sent to the back)
    */
  def dequeue: DequeueResult = {
    queue.dequeueOption match {
      case Some(actorAndNewQueue) =>
        val (placeholder, newQueue) = actorAndNewQueue
        val poolAcquisitionResult = pool.tryAcquire(placeholder.hogGroup)
        shuffleQueueAndGetLease(poolAcquisitionResult, queue, newQueue) match {
          case ShuffledQueueAndLease(lease, newQueue) => DequeueResult(newQueue, lease)
        }


        DequeueResult(Some(LeasedActor(placeholder.actor, shuffledQueueAndLease._1)), copy(queue = newQueue))
    }
  }

  def shuffleQueueAndGetLease(result: UnhoggableTokenPoolResult,
                              oldQueue: Queue[TokenQueuePlaceholder],
                              newQueue: Queue[TokenQueuePlaceholder]): ShuffledQueueAndLease = {
???
  }

  /**
    * Returns true if there's at least on element that can be dequeued, false otherwise.
    */
  def available: Boolean = queue.nonEmpty && pool.leased() < pool.capacity
}
