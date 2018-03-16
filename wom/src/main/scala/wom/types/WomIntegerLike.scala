package wom.types


object WomIntegerLike {
  implicit class EnhancedLong(val long: Long) extends AnyVal {
    def inIntRange: Boolean = long >= Int.MinValue && long <= Int.MaxValue
  }
}
