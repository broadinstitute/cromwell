package common.util

import cats.effect.IO

object IOUtil {

  def toIO[A](option: Option[A], errorMsg: String): IO[A] = {
    IO.fromEither(option.toRight(new RuntimeException(errorMsg)))
  }
}
