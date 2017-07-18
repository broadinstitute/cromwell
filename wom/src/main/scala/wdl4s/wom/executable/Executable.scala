package wdl4s.wom.executable

import wdl4s.wom.callable.Callable

/**
  * Closely related to the WdlNamespace, contains a set of Workflows and Tasks with a single Callable selected as the
  * entry point.
  */
final case class Executable(entryPoint: Callable) {
  val graph = entryPoint.graph
}
