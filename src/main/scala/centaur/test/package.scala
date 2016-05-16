package centaur
import cats.data.ValidatedNel

package object test {
  type ErrorOr[+A] = ValidatedNel[String, A]
}
