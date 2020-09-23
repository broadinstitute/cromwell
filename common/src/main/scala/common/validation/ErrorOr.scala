package common.validation

import cats.data.Validated.{Invalid, Valid}
import cats.data.{NonEmptyList, Validated}
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.instances.list._

object ErrorOr {
  type ErrorOr[+A] = Validated[NonEmptyList[String], A]

  implicit class EnhancedErrorOr[A](val eoa: ErrorOr[A]) extends AnyVal {
    def contextualizeErrors(s: => String): ErrorOr[A] = eoa.leftMap { errors =>
      val total = errors.size
      errors.zipWithIndex map { case (e, i) => s"Failed to $s (reason ${i + 1} of $total): $e" }
    }
  }

  implicit class ShortCircuitingFlatMap[A](val fa: ErrorOr[A]) extends AnyVal {
    /**
      * Not consistent with `Applicative#ap` but useful in for comprehensions.
      *
      * @see http://typelevel.org/cats/tut/validated.html#of-flatmaps-and-xors
      */
    def flatMap[B](f: A => ErrorOr[B]): ErrorOr[B] = {
      fa match {
        case Valid(v) => f(v)
        case i @ Invalid(_) => i
      }
    }
  }

  implicit class MapErrorOrRhs[A, B](val m: Map[A, ErrorOr[B]]) extends AnyVal {
    def sequence: ErrorOr[Map[A, B]] = m.traverseValues(identity)
  }

  implicit class MapTraversal[A, B](val m: Map[A, B]) extends AnyVal {
    def traverse[C,D](f: ((A,B)) => ErrorOr[(C,D)]): ErrorOr[Map[C,D]] = m.toList.traverse(f).map(_.toMap)
    def traverseValues[C](f: B => ErrorOr[C]): ErrorOr[Map[A,C]] = m.traverse { case (a, b) => f(b).map(c => (a,c)) }
  }

  // Note! See the bottom of this file for a generator function for 2 through 22 of these near-identical ShortCircuitingFlatMapTupleNs...

  implicit class ShortCircuitingFlatMapTuple1[A](val t1: Tuple1[ErrorOr[A]]) extends AnyVal {
    def flatMapN[T_OUT](f_1: (A) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t1._1.flatMap(f_1)
  }

  implicit class ShortCircuitingFlatMapTuple2[A, B](val t2: (ErrorOr[A], ErrorOr[B])) extends AnyVal {
    def flatMapN[T_OUT](f2: (A, B) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t2.tupled flatMap f2.tupled
  }

  implicit class ShortCircuitingFlatMapTuple3[A, B, C](val t3: (ErrorOr[A], ErrorOr[B], ErrorOr[C])) extends AnyVal {
    def flatMapN[T_OUT](f3: (A, B, C) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t3.tupled flatMap f3.tupled
  }

  implicit class ShortCircuitingFlatMapTuple4[A, B, C, D](val t4: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D])) extends AnyVal {
    def flatMapN[T_OUT](f4: (A, B, C, D) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t4.tupled flatMap f4.tupled
  }

  implicit class ShortCircuitingFlatMapTuple5[A, B, C, D, E](val t5: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E])) extends AnyVal {
    def flatMapN[T_OUT](f5: (A, B, C, D, E) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t5.tupled flatMap f5.tupled
  }

  implicit class ShortCircuitingFlatMapTuple6[A, B, C, D, E, F](val t6: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F])) extends AnyVal {
    def flatMapN[T_OUT](f6: (A, B, C, D, E, F) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t6.tupled flatMap f6.tupled
  }

  implicit class ShortCircuitingFlatMapTuple7[A, B, C, D, E, F, G](val t7: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G])) extends AnyVal {
    def flatMapN[T_OUT](f7: (A, B, C, D, E, F, G) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t7.tupled flatMap f7.tupled
  }

  implicit class ShortCircuitingFlatMapTuple8[A, B, C, D, E, F, G, H](val t8: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H])) extends AnyVal {
    def flatMapN[T_OUT](f8: (A, B, C, D, E, F, G, H) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t8.tupled flatMap f8.tupled
  }

  implicit class ShortCircuitingFlatMapTuple9[A, B, C, D, E, F, G, H, I](val t9: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I])) extends AnyVal {
    def flatMapN[T_OUT](f9: (A, B, C, D, E, F, G, H, I) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t9.tupled flatMap f9.tupled
  }

  implicit class ShortCircuitingFlatMapTuple10[A, B, C, D, E, F, G, H, I, J](val t10: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J])) extends AnyVal {
    def flatMapN[T_OUT](f10: (A, B, C, D, E, F, G, H, I, J) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t10.tupled flatMap f10.tupled
  }

  implicit class ShortCircuitingFlatMapTuple11[A, B, C, D, E, F, G, H, I, J, K](val t11: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K])) extends AnyVal {
    def flatMapN[T_OUT](f11: (A, B, C, D, E, F, G, H, I, J, K) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t11.tupled flatMap f11.tupled
  }

  implicit class ShortCircuitingFlatMapTuple12[A, B, C, D, E, F, G, H, I, J, K, L](val t12: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L])) extends AnyVal {
    def flatMapN[T_OUT](f12: (A, B, C, D, E, F, G, H, I, J, K, L) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t12.tupled flatMap f12.tupled
  }

  implicit class ShortCircuitingFlatMapTuple13[A, B, C, D, E, F, G, H, I, J, K, L, M](val t13: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M])) extends AnyVal {
    def flatMapN[T_OUT](f13: (A, B, C, D, E, F, G, H, I, J, K, L, M) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t13.tupled flatMap f13.tupled
  }

  implicit class ShortCircuitingFlatMapTuple14[A, B, C, D, E, F, G, H, I, J, K, L, M, N](val t14: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N])) extends AnyVal {
    def flatMapN[T_OUT](f14: (A, B, C, D, E, F, G, H, I, J, K, L, M, N) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t14.tupled flatMap f14.tupled
  }

  implicit class ShortCircuitingFlatMapTuple15[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O](val t15: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O])) extends AnyVal {
    def flatMapN[T_OUT](f15: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t15.tupled flatMap f15.tupled
  }

  implicit class ShortCircuitingFlatMapTuple16[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P](val t16: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P])) extends AnyVal {
    def flatMapN[T_OUT](f16: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t16.tupled flatMap f16.tupled
  }

  implicit class ShortCircuitingFlatMapTuple17[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q](val t17: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P], ErrorOr[Q])) extends AnyVal {
    def flatMapN[T_OUT](f17: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t17.tupled flatMap f17.tupled
  }

  implicit class ShortCircuitingFlatMapTuple18[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R](val t18: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P], ErrorOr[Q], ErrorOr[R])) extends AnyVal {
    def flatMapN[T_OUT](f18: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t18.tupled flatMap f18.tupled
  }

  implicit class ShortCircuitingFlatMapTuple19[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S](val t19: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P], ErrorOr[Q], ErrorOr[R], ErrorOr[S])) extends AnyVal {
    def flatMapN[T_OUT](f19: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t19.tupled flatMap f19.tupled
  }

  implicit class ShortCircuitingFlatMapTuple20[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T](val t20: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P], ErrorOr[Q], ErrorOr[R], ErrorOr[S], ErrorOr[T])) extends AnyVal {
    def flatMapN[T_OUT](f20: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t20.tupled flatMap f20.tupled
  }

  implicit class ShortCircuitingFlatMapTuple21[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U](val t21: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P], ErrorOr[Q], ErrorOr[R], ErrorOr[S], ErrorOr[T], ErrorOr[U])) extends AnyVal {
    def flatMapN[T_OUT](f21: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t21.tupled flatMap f21.tupled
  }

  implicit class ShortCircuitingFlatMapTuple22[A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V](val t22: (ErrorOr[A], ErrorOr[B], ErrorOr[C], ErrorOr[D], ErrorOr[E], ErrorOr[F], ErrorOr[G], ErrorOr[H], ErrorOr[I], ErrorOr[J], ErrorOr[K], ErrorOr[L], ErrorOr[M], ErrorOr[N], ErrorOr[O], ErrorOr[P], ErrorOr[Q], ErrorOr[R], ErrorOr[S], ErrorOr[T], ErrorOr[U], ErrorOr[V])) extends AnyVal {
    def flatMapN[T_OUT](f22: (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t22.tupled flatMap f22.tupled
  }

  object ErrorOrGen {
    /**
      * Because maintaining 22 near-identical functions is a bore...
      * This function can regenerate them if we need to make changes.
      * Example usage:
      *   System.out.println((2 to 22).map(mkShortCircuitingFlatMapTupleNFunction).mkString(""))
      */
    def mkShortCircuitingFlatMapTupleNFunction(n: Int): String = {

      def getCharForNumber(i: Int): String = String.valueOf((i + 65).asInstanceOf[Char])
      def abcList(n: Int): Seq[String] = (0 until n).map(getCharForNumber)

      def aList(n: Int): String = abcList(n).mkString(", ")
      def errorOrList(n: Int): String = abcList(n).map(a => s"ErrorOr[$a]").mkString(", ")

      def line1(n: Int): String =
        s"implicit class ShortCircuitingFlatMapTuple$n[${aList(n)}](val t$n: (${errorOrList(n)})) extends AnyVal {"

      def line2(n: Int): String =
        s"  def flatMapN[T_OUT](f$n: (${aList(n)}) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t$n.tupled flatMap f$n.tupled"

      s"""${line1(n)}
         |  ${line2(n)}
         |}
         |
       |""".stripMargin
    }
  }
}

