package wdl4s.wom.executable

import lenthall.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl4s.wdl.WorkflowRawInputs
import wdl4s.wom.callable.Callable
import wdl4s.wom.graph.Graph.ResolvedWorkflowInput
import wdl4s.wom.graph.GraphNodePort.OutputPort

/**
  * Closely related to the WdlNamespace, contains a set of Workflows and Tasks with a single Callable selected as the
  * entry point.
  */
final case class Executable(entryPoint: Callable) {
  val graph = entryPoint.graph
  
  def validateWorkflowInputs(inputsMapping: WorkflowRawInputs): ErrorOr[Map[OutputPort, ResolvedWorkflowInput]] = for {
    validGraph <- graph
    inputs <- validGraph.validateWorkflowInputs(inputsMapping)
  } yield inputs
}
