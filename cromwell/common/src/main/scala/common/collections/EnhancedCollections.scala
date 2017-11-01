package common.collections

import scala.collection.TraversableLike
import scala.collection.generic.CanBuildFrom
import scala.reflect.ClassTag

object EnhancedCollections {

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
}
