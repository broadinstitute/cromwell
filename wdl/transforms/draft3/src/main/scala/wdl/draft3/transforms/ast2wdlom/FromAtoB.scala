package wdl.draft3.transforms.ast2wdlom

import cats.syntax.traverse._
import cats.instances.vector._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr.ShortCircuitingFlatMap

trait FromAtoB[A, B] {
  def convert(a: A): ErrorOr[B]
}

object FromAtoB {
  implicit def apply[A, B](implicit fromImplicitScope: FromAtoB[A, B]) = fromImplicitScope

  def viaX[A, X, B](implicit ax: FromAtoB[A, X], xb: FromAtoB[X, B]): FromAtoB[A, B] = new FromAtoB[A, B] {
    override def convert(a: A): ErrorOr[B] = { ax.convert(a) flatMap xb.convert }
  }
  def forVectors[A, B](implicit ab: FromAtoB[A, B]): FromAtoB[Vector[A], Vector[B]] = new FromAtoB[Vector[A], Vector[B]] {
    override def convert(a: Vector[A]): ErrorOr[Vector[B]] = a.traverse[ErrorOr, B] { ab.convert }
  }
}
