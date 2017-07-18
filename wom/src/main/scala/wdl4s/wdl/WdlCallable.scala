package wdl4s.wdl

import wdl4s.wom.callable.Callable

trait WdlCallable extends Scope {
  def womDefinition: Callable
  def outputs: Seq[Output]
}
