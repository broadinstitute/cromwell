package cromwell.engine.workflow.tokens

import akka.actor.ActorRef
import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.TokenQueue.{DequeuedActor, LeasedActor}
import io.github.andrebeat.pool.Lease

import scala.collection.immutable.Queue

object TokenQueue {
  case class DequeuedActor(leasedActor: LeasedActor, tokenQueue: TokenQueue)
  case class LeasedActor(actor: ActorRef, lease: Lease[JobExecutionToken])
  def apply(tokenType: JobExecutionTokenType) = new TokenQueue(Queue.empty, TokenPool(tokenType))
}

/**
  * A queue assigned to a pool.
  * Elements can be dequeued if the queue is not empty and there are tokens in the pool
  */
final case class TokenQueue(queue: Queue[ActorRef], private [tokens] val pool: TokenPool) {
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
  def enqueue(actor: ActorRef): TokenQueue = {
    copy(queue = queue.enqueue(actor))
  }

  /**
    * Returns a DequeuedActor if there's an element available, None otherwise.
    */
  def dequeueOption: Option[DequeuedActor] = {
    for {
      actorAndNewQueue <- queue.dequeueOption
      (actor, newQueue) = actorAndNewQueue
      lease <- pool.tryAcquire()
    } yield DequeuedActor(LeasedActor(actor, lease), copy(queue = newQueue))
  }

  /**
    * Returns true if there's at least on element that can be dequeued, false otherwise.
    */
  def available: Boolean = queue.nonEmpty && pool.leased() < pool.capacity
}
