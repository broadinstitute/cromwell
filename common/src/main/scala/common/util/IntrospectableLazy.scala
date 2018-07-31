package common.util

import scala.language.implicitConversions

// Adapted from https://stackoverflow.com/a/17057490/818054

object IntrospectableLazy {

  def lazily[A](f: => A): IntrospectableLazy[A] = new IntrospectableLazy(f)

  implicit def evalLazy[A](l: IntrospectableLazy[A]): A = l()

}

class IntrospectableLazy[A] private(f: => A) {

  private var option: Option[A] = None

  def apply(): A = {
    option match {
      case Some(a) => a
      case None =>
        synchronized { option match {
          case Some(a) => a
          case None => val a = f; option = Some(a); a
        }}
    }
  }

  def exists: Boolean = option.isDefined

}
