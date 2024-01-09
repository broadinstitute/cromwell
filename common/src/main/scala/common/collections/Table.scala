package common.collections

import scala.collection.immutable

object Table {

  /**
    * Instantiates an empty table
    */
  def empty[R, C, V]: Table[R, C, V] = Table(Map.empty[R, Map[C, V]])

  /**
    * Fills a new table with the values
    */
  def fill[R, C, V](values: Iterable[(R, C, V)]): Table[R, C, V] = empty[R, C, V].addAll(values)
}

/**
  * Data structure to represent a table (with rows, columns and values)
  * Really just an alias for Map[R, Map[C, V]]
  * @tparam R type of the row
  * @tparam C type of the column
  * @tparam V type of the value
  */
case class Table[R, C, V](table: Map[R, Map[C, V]]) {

  /**
    * Returns true if the table contains a value at row / column
    */
  def contains(row: R, column: C): Boolean = table.get(row).exists(_.contains(column))

  /**
    * Get the value at row / column
    */
  def getValue(row: R, column: C): Option[V] = table.get(row).flatMap(_.get(column))

  /**
    * Get the values at row if there is an entry, None otherwise
    */
  def rowOptional(row: R): Option[Map[C, V]] = table.get(row)

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
  def add(row: R, column: C, value: V): Table[R, C, V] =
    this.copy(
      table = table.updated(row, table.getOrElse(row, Map.empty).updated(column, value))
    )

  /**
    * Add all values
    */
  def addAll(values: Iterable[(R, C, V)]): Table[R, C, V] =
    values.foldLeft(this)(_.addTriplet(_))

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
