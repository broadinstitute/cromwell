package wom.types


object WomIntegerLike {
  implicit class EnhancedLong(val long: Long) extends AnyVal {
    def inIntRange: Boolean = long >= Int.MinValue && long <= Int.MaxValue
  }
  implicit class EnhancedDouble(val double: Double) extends AnyVal {
    def inIntRange: Boolean = double >= Int.MinValue && double <= Int.MaxValue
  }
}
