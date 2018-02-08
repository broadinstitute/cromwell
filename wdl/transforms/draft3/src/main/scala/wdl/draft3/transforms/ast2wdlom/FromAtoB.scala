package wdl.draft3.transforms.ast2wdlom

import common.validation.ErrorOr.ErrorOr

trait FromAtoB[A, B] {
  def convert(a: A): ErrorOr[B]
}
