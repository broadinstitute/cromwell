package cromwell.engine.workflow.lifecycle.execution.job.preparation

import akka.actor.Props
import common.validation.ErrorOr._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import wom.expression.IoFunctionSet
import wom.graph.CallNode
import wom.values.WomEvaluatedCallInputs

object CallPreparation {
  sealed trait CallPreparationActorCommands
  case class Start(valueStore: ValueStore) extends CallPreparationActorCommands

  trait CallPreparationActorResponse

  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props)
      extends CallPreparationActorResponse

  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse

  def resolveAndEvaluateInputs(callKey: CallKey,
                               expressionLanguageFunctions: IoFunctionSet,
                               valueStore: ValueStore
  ): ErrorOr[WomEvaluatedCallInputs] =
    CallNode.resolveAndEvaluateInputs(callKey.node, expressionLanguageFunctions, valueStore.resolve(callKey.index))
}
