import cats.data.NonEmptyList

package object lenthall {
  type Checked[A] = Either[NonEmptyList[String], A]

}
