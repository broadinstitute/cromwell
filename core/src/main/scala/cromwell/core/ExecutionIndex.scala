package cromwell.core

/*
 * Type and implicit conversion classes for ExecutionIndex
 */
object ExecutionIndex {
  type ExecutionIndex = Option[Int]
  val IndexNone = -1 // "It's a feature" https://bugs.mysql.com/bug.php?id=8173

  implicit class IndexEnhancedInt(val value: Int) extends AnyVal {
    def toIndex: ExecutionIndex = value match {
      case IndexNone => None
      case i => Option(i)
    }
  }

  implicit class IndexEnhancedIndex(val index: ExecutionIndex) extends AnyVal {
    def fromIndex: Int = index match {
      case None => IndexNone
      case Some(i) => i
    }
    def isShard: Boolean = index.nonEmpty
  }

  implicit val ExecutionIndexOrdering = new Ordering[ExecutionIndex] {
    override def compare(x: ExecutionIndex, y: ExecutionIndex): Int = {
      x.fromIndex.compareTo(y.fromIndex)
    }
  }
}
