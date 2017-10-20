package wdl

import common.validation.ErrorOr.ErrorOr
import wom.callable.Callable

trait WdlCallable extends Scope {
  def womDefinition: ErrorOr[Callable]
  def outputs: Seq[Output]
}
