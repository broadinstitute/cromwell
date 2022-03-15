package common

package object assertion {
  object ManyTimes {
    implicit class intWithTimes(n: Int) {
      def times(f: => Unit) = 1 to n foreach { _ => f }
      def indexedTimes(f: Int => Unit) = 0 until n foreach { i => f(i) }
      def of[A](f: () => A) = (0 until n map { _ => f() }).toVector
    }
  }
}
