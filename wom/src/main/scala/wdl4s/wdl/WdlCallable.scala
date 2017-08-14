package wdl4s.wdl

import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.callable.Callable

trait WdlCallable extends Scope {
  def womDefinition: ErrorOr[Callable]
  def outputs: Seq[Output]
}
