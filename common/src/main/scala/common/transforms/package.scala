package common

import cats.data.Kleisli

import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import common.validation.Checked._

package object transforms {
  type CheckedAtoB[A, B] = Kleisli[Checked, A, B]

  object CheckedAtoB {
    def apply[A, B](implicit runner: CheckedAtoB[A, B]): CheckedAtoB[A, B] = runner
    def fromCheck[A, B](run: A => Checked[B]): CheckedAtoB[A, B] = Kleisli(run)
    def fromCheck[A, B](context: String)(run: A => Checked[B]): CheckedAtoB[A, B] = Kleisli(
      runCheckWithContext(run, _ => context)
    )
    def fromCheck[A, B](context: A => String)(run: A => Checked[B]): CheckedAtoB[A, B] = Kleisli(
      runCheckWithContext(run, context)
    )
    def fromErrorOr[A, B](run: A => ErrorOr[B]): CheckedAtoB[A, B] = Kleisli(runThenCheck(run))
    def fromErrorOr[A, B](context: String)(run: A => ErrorOr[B]): CheckedAtoB[A, B] = Kleisli(
      runErrorOrWithContext(run, _ => context)
    )
    def fromErrorOr[A, B](context: A => String)(run: A => ErrorOr[B]): CheckedAtoB[A, B] = Kleisli(
      runErrorOrWithContext(run, context)
    )
    private def runThenCheck[A, B](run: A => ErrorOr[B]): A => Checked[B] = (a: A) => run(a).toEither
    private def runErrorOrWithContext[A, B](run: A => ErrorOr[B], context: A => String): A => Checked[B] = (a: A) =>
      run(a).toEither.contextualizeErrors(context(a))
    private def runCheckWithContext[A, B](run: A => Checked[B], context: A => String): A => Checked[B] = (a: A) =>
      run(a).contextualizeErrors(context(a))

    def firstSuccess[A, B](options: List[CheckedAtoB[A, B]], operationName: String): CheckedAtoB[A, B] =
      Kleisli[Checked, A, B] { a =>
        if (options.isEmpty) {
          s"Unable to $operationName: No import resolvers provided".invalidNelCheck
        } else {
          val firstAttempt = options.head.run(a)
          options.tail.foldLeft[Checked[B]](firstAttempt) { (currentResult, nextOption) =>
            currentResult match {
              case v: Right[_, _] => v
              case Left(currentErrors) =>
                nextOption.run(a) match {
                  case v: Right[_, _] => v
                  case Left(newErrors) => Left(currentErrors ++ newErrors.toList)
                }
            }
          }
        }
      }
  }
}
