package common.collections

import common.collections.EnhancedCollections._
import scala.collection.immutable.Queue

object WeightedQueue {
  def empty[T, W](weightFunction: T => W)(implicit n: Numeric[W]) = {
    WeightedQueue(Queue.empty[T], weightFunction, n.zero)
  }
}

/**
  * Wrapper class around an immutable queue.
  * In addition to the queue, a weight function is provided that provides the weight W of an element T.
  * The total weight of the queue is accessible, as well as a method to take the head of the queue based on a max weight value.
  */
final case class WeightedQueue[T, W](innerQueue: Queue[T],
                                     private val weightFunction: T => W,
                                     weight: W)(implicit n: Numeric[W]) {
  import n._

  def enqueue(element: T): WeightedQueue[T, W] = {
    this.copy(innerQueue = innerQueue.enqueue(element), weight = weight + weightFunction(element))
  }

  def dequeue: (T, WeightedQueue[T, W]) = {
    val (element, tail) = innerQueue.dequeue
    element -> this.copy(innerQueue = tail, weight = weight - weightFunction(element))
  }

  def dequeueOption: Option[(T, WeightedQueue[T, W])] = {
    innerQueue.dequeueOption map {
      case (element, tail) =>
        element -> this.copy(innerQueue = tail, weight = weight - weightFunction(element))
    }
  }
  
  def behead(maxWeight: W, maxLength: Option[Int] = None, strict: Boolean = false): (Vector[T], WeightedQueue[T, W]) = {
    val DeQueued(head, tail) = innerQueue.takeWhileWeighted(maxWeight, weightFunction, maxLength, strict)
    head -> this.copy(innerQueue = tail, weight = weight - head.map(weightFunction).sum)
  }
}
