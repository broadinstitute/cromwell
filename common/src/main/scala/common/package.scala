import cats.data.NonEmptyList

package object common {
  type Checked[+A] = Either[NonEmptyList[String], A]
}
