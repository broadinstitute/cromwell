package cromwell.engine.workflow.tokens

import cromwell.engine.workflow.tokens.TokenQueue.{DequeueResult, LeasedActor}

/**
  * Creates an Iterator from a list of TokenQueues.
  * It will keep rotating the list until it finds a queue with an element that can be dequeued.
  * If no queue can be dequeued, the iterator is empty.
  */
final class RoundRobinQueueIterator(initialTokenQueue: List[TokenQueue], initialPointer: Int) extends Iterator[LeasedActor] {
  // Assumes the number of queues won't change during iteration (it really shouldn't !)
  private val numberOfQueues = initialTokenQueue.size
  // Indicate the index of next queue to try to dequeue from.
  // In case token queues have been removed since the last time an iterator was created on this token queue list, make sure the pointer is in range.
  private[this] var pointer: Int = if (numberOfQueues == 0) 0 else initialPointer % numberOfQueues
  // List of queues available
  private[this] var tokenQueues: List[TokenQueue] = initialTokenQueue

  /**
    * As we iterate and dequeue elements, and because the queues are immutable,
    * they need to be kept up to date with their new version every time an element is dequeued.
    * This returns the up to date version of the queues according to the current state of iteration.
    */
  def updatedQueues = tokenQueues

  /**
    * And the same w.r.t the current pointer:
    */
  def updatedPointer = pointer

  def hasNext = tokenQueues.exists(_.available)
  def next() = findFirst.getOrElse(unexpectedlyEmpty)

  def unexpectedlyEmpty: LeasedActor =
    throw new IllegalStateException("Token iterator is empty")

  // Goes over the JobExecutionTokenType queues and returns the first element that can be dequeued.
  // NB: In normal usage, there's one queue per backend.
  private def findFirst: Option[LeasedActor] = {
    // Starting from pointer make a list of indices circling back to pointer
    // This will ensure we try all the queues, as well as keep rotating
    // and don't keep emptying the same queue as long as it has elements
    // For instance, if we have 5 queues and pointer is 2, we want to try indices (2, 3, 4, 0, 1)

    val indexStream = ((pointer until numberOfQueues) ++ (0 until pointer)).toStream
    val dequeuedTokenStream = indexStream.map(index => tokenQueues(index).dequeue -> index)

    val firstLeasedActor = dequeuedTokenStream.collectFirst({
      case (DequeueResult(Some(dequeuedActor), newTokenQueue), index) =>
        // Update the tokenQueues with the new queue
        tokenQueues = tokenQueues.updated(index, newTokenQueue)
        // Update the index. Add 1 to force trying all the queues as we call next, even if the first one is available
        pointer = (index + 1) % numberOfQueues
        dequeuedActor
    })

    firstLeasedActor
  }
}
