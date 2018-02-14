package common

import cats.data.Kleisli
import common.validation.ErrorOr.ErrorOr

package object transforms {
  type CheckedAtoB[A, B] = Kleisli[Checked, A, B]

  object CheckedAtoB {
    def apply[A, B](implicit runner: CheckedAtoB[A, B]): CheckedAtoB[A, B] = runner
    def fromCheck[A, B](run: A => Checked[B]): CheckedAtoB[A, B] = Kleisli(run)
    def fromErrorOr[A, B](run: A => ErrorOr[B]): CheckedAtoB[A, B] = Kleisli(runThenCheck(run))

    private def runThenCheck[A, B](run: A => ErrorOr[B]): A => Checked[B] = (a: A) => { run(a).toEither }
  }
}
