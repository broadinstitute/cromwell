package lenthall.collections

import lenthall.collections.Table.Table

import scala.collection.{TraversableLike, immutable}
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
  
  implicit class EnhancedTable[R, C, V](val table: Table[R, C, V]) extends AnyVal {
    /**
      * Returns true if the table contains a value at row / column
      */
    def contains(row: R, column: C): Boolean = table.get(row).exists(_.contains(column))

    /**
      * Get the value at row / column
      */
    def get(row: R, column: C): Option[V] = table.get(row).flatMap(_.get(column))

    /**
      * Get the row
      */
    def row(row: R): Map[C, V] = table.getOrElse(row, Map.empty)

    /**
      * Get the column (more expansive than row)
      */
    def column(column: C): Map[R, V] = table flatMap { case (rowKey, col) => col.get(column).map(rowKey -> _) }

    /**
      * Add a value at row / column
      */
    def add(row: R, column: C, value: V): Table[R, C, V] = table.updated(row, table.getOrElse(row, Map.empty).updated(column, value))

    /**
      * Add all values
      */
    def addAll(values: Iterable[(R, C, V)]): Table[R, C, V] = values.foldLeft(table)(_.addTriplet(_))

    /**
      * Add a value as a triplet
      */
    def addTriplet(value: (R, C, V)): Table[R, C, V] = Function.tupled(add _)(value)

    /**
      * Returns all values as triplets
      */
    def valuesTriplet: immutable.Iterable[(R, C, V)] = for {
      (row, columns) <- table
      (column, value) <- columns
    } yield (row, column, value)
  }
}

