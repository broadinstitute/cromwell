package cromwell.engine.workflow.tokens

import java.util.concurrent.atomic.AtomicBoolean

import akka.actor.ActorRef
import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.TokenQueue.{DequeuedActor, LeasedActor, TokenQueuePlaceholder}
import io.github.andrebeat.pool.Lease

import scala.collection.immutable.Queue

object TokenQueue {
  final case class DequeuedActor(leasedActor: LeasedActor, tokenQueue: TokenQueue)
  final case class LeasedActor(actor: ActorRef, lease: Lease[JobExecutionToken])
  def apply(tokenType: JobExecutionTokenType) = new TokenQueue(Queue.empty, TokenPool(tokenType))
  final case class TokenQueuePlaceholder(actor: ActorRef, hogGroup: String)
  final case class TokenQueueState(queue: Queue[TokenQueuePlaceholder], hogCounts: Map[String, Int])

  final case class TokenHoggingLease(lease: Lease[JobExecutionToken]) extends Lease[JobExecutionToken] {
    private[this] val dirty = new AtomicBoolean(false)
    override protected[this] def a: JobExecutionToken = lease.get

    override protected[this] def handleRelease(): Unit = {
      if (dirty.compareAndSet(expect = false, update = true)) {

      }
      lease.release()
    }
    override protected[this] def handleInvalidate(): Unit = lease.invalidate()
  }
}

/**
  * A queue assigned to a pool.
  * Elements can be dequeued if the queue is not empty and there are tokens in the pool
  */
final case class TokenQueue(queue: Queue[TokenQueuePlaceholder], hogIndex: Int, private [tokens] val pool: TokenPool) {
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
    * Returns Some(DequeuedActor(...)) if there's:
    *  - an actor to dequeue
    *  - a token available for it
    *  Returns None otherwise.
    */
  def dequeueOption: Option[DequeuedActor] = {
    for {
      actorAndNewQueue <- queue.dequeueOption
      (queuePlaceholder, newQueue) = actorAndNewQueue
      lease <- acquireLease(queuePlaceholder.hogGroup, pool)
    } yield DequeuedActor(LeasedActor(queuePlaceholder.actor, lease), copy(queue = newQueue))
  }

  def acquireLease(hogGroup: String, pool: TokenPool): Option[Lease[JobExecutionToken]] = {
    pool.tryAcquire()
  }

  /**
    * Returns true if there's at least on element that can be dequeued, false otherwise.
    */
  def available: Boolean = queue.nonEmpty && pool.leased() < pool.capacity
}

