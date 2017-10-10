package lenthall.collections

object Table {
  import lenthall.collections.EnhancedCollections._

  /**
    * Data structure to represent a table (with rows, columns and values)
    * Really just an alias for Map[R, Map[C, V]]
    * @tparam R type of the row
    * @tparam C type of the column
    * @tparam V type of the value
    */
  type Table[R, C, V] = Map[R, Map[C, V]]

  /**
    * Instantiates an empty table
    */
  def empty[R, C, V]: Table[R, C, V] = Map.empty[R, Map[C, V]]

  /**
    * Fills a new table with the values
    */
  def fill[R, C, V](values: Iterable[(R, C, V)]): Table[R, C, V] = empty[R, C, V].addAll(values)
}
