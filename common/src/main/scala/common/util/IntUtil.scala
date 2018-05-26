package common.util

object IntUtil {
  implicit class EnhancedInt(val int: Int) extends AnyVal {
    def isEven: Boolean = int % 2 == 0
    def isOdd: Boolean = !isEven
  }
}
