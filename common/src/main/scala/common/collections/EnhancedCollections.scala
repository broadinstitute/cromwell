package common.collections

import scala.annotation.tailrec
import scala.collection.TraversableLike
import scala.collection.generic.CanBuildFrom
import scala.collection.immutable.{MapLike, Queue}
import scala.reflect.ClassTag

object EnhancedCollections {
  
  case class DeQueued[A](head: Vector[A], tail: Queue[A])

  /**
    * After trying and failing to do this myself, I got this to work by copying the answer from here:
    * https://stackoverflow.com/questions/29886246/scala-filter-by-type
    */
  implicit class EnhancedTraversableLike[T2, Repr <: TraversableLike[T2, Repr], That](val traversable: TraversableLike[T2, Repr]) extends AnyVal {
    /**
      * Lets you filter a collection by type.
      *
      * Warning: intelliJ has problems working out the return type but it is what you'd expect it to be.
      * If you dislike intelliJ red, you can use type ascription to give it a hand
      *
      * eg.
      * val xs: Set[Object]
      * val strings: Set[String] = xs.filterByType[String]
      */
    def filterByType[T <: T2](implicit tag: ClassTag[T], bf: CanBuildFrom[Repr, T, That]): That = traversable.collect { case t: T => t }

    def firstByType[T <: T2](implicit tag: ClassTag[T]): Option[T] = traversable collectFirst { case t: T => t }
  }

  implicit class EnhancedQueue[A](val queue: Queue[A]) extends AnyVal {

    /**
      * Take from the queue based on the following conditions:
      * (in this context, "head" means the Vector[A] dequeued, not a single A element)
      * 
      * The head's total weight (calculated using weightFunction) will be:
      *   - Always less than or equal to maxWeight if strict is true
      *   - The smallest list of elements at the head of the queue for which the total weight is greater or equal to max weight if strict is false
      *   
      * In other words, if a single element happens to be greater than max weight:
      *   with strict = true the element will be dropped and processing will continue
      *   with strict = false the element will be added to the head and processing will stop
      * 
      * If maxHeadLength is provided, the head size will be <= maxHeadLength's value.
      * 
      * --- Implementation ---
      * The implementation is simple but naive. The more general problem is described here https://en.wikipedia.org/wiki/Knapsack_problem
      * 
      * @tparam W a numeric type for the weight
      * @return a tuple of 
      *         The elements removed from the queue as a vector
      *         The rest of the queue (tail)
      */
    def takeWhileWeighted[W](maxWeight: W,
                             weightFunction: A => W,
                             maxHeadLength: Option[Int],
                             strict: Boolean = false)
                            (implicit n: Numeric[W], c: Ordering[W]): DeQueued[A] = {
      import n._

      @tailrec
      def takeWhileWeightedRec(tail: Queue[A], head: Vector[A], weight: W): (Vector[A], Queue[A]) = {
        // Stay under maxHeadLength if it's specified
        if (maxHeadLength.exists(head.size >= _)) head -> tail
        else tail.dequeueOption
          .map({
            // Compute the dequeued element's weight
            case (element, dequeued) => (element, weightFunction(element), dequeued)
          }) match {
          // If the element's weight is > maxWeight and strict is true, drop the element
          case Some((_, elementWeight, dequeued)) if c.gteq(elementWeight, maxWeight) && strict => takeWhileWeightedRec(dequeued, head, weight)
          // If we're under the max weight, add the element to the head and recurse
          case Some((element, elementWeight, dequeued)) if c.lteq(elementWeight + weight, maxWeight) => takeWhileWeightedRec(dequeued, head :+ element, weight + elementWeight)
          // Otherwise stop here (make sure to return the original queue so we don't lose the last dequeued element)
          case _ => head -> tail
        }
      }

      if (queue.isEmpty || maxHeadLength.contains(0)) DeQueued(Vector.empty, queue)
      // If strict is enabled, we should never return a head with a weight > maxWeight. So start from the original queue and drop elements over maxWeight if necessary 
      else if (strict) {
        val (head, tail) = takeWhileWeightedRec(queue, Vector.empty, n.zero)
        DeQueued(head, tail)
      }
      // Otherwise to ensure we don't deadlock, start the recursion with the head of the queue, this way even if it's over maxWeight it'll return a single element head 
      else {
        val (head, tail) = takeWhileWeightedRec(queue.tail, queue.headOption.toVector, queue.headOption.map(weightFunction).getOrElse(n.zero))
        DeQueued(head, tail)
      }
    }
  }

  implicit class EnhancedMapLike[A, +B, +This <: MapLike[A, B, This] with Map[A, B]](val mapLike: MapLike[A, B, This]) {
    /**
      * 'safe' in that unlike the implementation hiding behind `MapLike#mapValues` this is strict. i.e. this will only
      * evaluate the supplied function once on each value and at the time this method is called.
      */
    def safeMapValues[C](f: B => C): Map[A, C] = mapLike map { case (k, v) => k -> f(v) }
  }
}
