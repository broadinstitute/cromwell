package common.legacy

import scala.util.{Failure, Right, Success, Try}

/**
  * Some stuff to support 2.11 that can go away if/when we only need to support 2.12:
  */
object TwoElevenSupport {

  /**
    * Emulates the 'toEither' method on Try in 2.12:
    */
  implicit class TacticallyEnhancedTry[A](t: Try[A]) {
    def tacticalToEither: Either[Throwable, A] = t match {
      case Success(a) => Right(a)
      case Failure(e) => Left(e)
    }
  }
}
