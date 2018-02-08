package wdl

import common.validation.ErrorOr.ErrorOr
import wom.callable.Callable

trait WdlCallable extends Scope {
  def outputs: Seq[Output]
  lazy val inputNames: Seq[String] = this match {
    case _: WdlWorkflow => declarations.filter(_.upstream.isEmpty).map(_.unqualifiedName)
    case _: WdlTask => declarations.map(_.unqualifiedName)
  }
}
