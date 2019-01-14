package common

import cats.effect.IO

package object exception {

  def toIO[A](option: Option[A], errorMsg: String): IO[A] = {
    IO.fromEither(option.toRight(new RuntimeException(errorMsg)))
  }
}
