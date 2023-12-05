package common.util

import scala.language.implicitConversions

// Designed to hold resources that may or may not be initialized during the lifetime
// of the application, much like lazy.
//
// The difference is that e.g. a shutdown routine can check whether an IntrospectableLazy
// actually exists() and skip cleaning it up if it doesn't.
//
// Otherwise, the shutdown may actually be the first reference to the resource, causing
// it to be initialized immediately before cleaning it up.

// Adapted from https://stackoverflow.com/a/17057490/818054

object IntrospectableLazy {

  def lazily[A](f: => A): IntrospectableLazy[A] = new IntrospectableLazy(f)

  implicit def evalLazy[A](l: IntrospectableLazy[A]): A = l()

}

class IntrospectableLazy[A] private (f: => A) {

  private var option: Option[A] = None

  def apply(): A =
    option match {
      case Some(a) => a
      case None =>
        synchronized {
          option match {
            case Some(a) => a
            case None => val a = f; option = Some(a); a
          }
        }
    }

  def exists: Boolean = option.isDefined

}
