package cromwell.engine.workflow.lifecycle.execution.keys

import akka.actor.ActorRef
import cats.syntax.validated._
import common.validation.ErrorOr._
import cromwell.core.ExecutionIndex._
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionDiff
import cromwell.engine.workflow.lifecycle.execution.keys.ExpressionKey.{ExpressionEvaluationFailedResponse, ExpressionEvaluationSucceededResponse}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExpressionNodeLike
import wom.values.WomValue

final case class ExpressionKey(node: ExpressionNodeLike, index: ExecutionIndex) extends JobKey {
  override val attempt = 1
  override lazy val tag = s"Expression-${node.localName}:${index.fromIndex}:$attempt"

  def processRunnable(ioFunctionSet: IoFunctionSet, valueStore: ValueStore, workflowExecutionActor: ActorRef): ErrorOr[WorkflowExecutionDiff] = {
    // Send a message to self in case we decide to change evaluate to return asynchronously, if we don't we could
    // directly add the value to the value store in the execution diff
    node.evaluate(valueStore.resolve(index), ioFunctionSet) match {
      case Right(result) => workflowExecutionActor ! ExpressionEvaluationSucceededResponse(this, result)
      case Left(f) => 
        workflowExecutionActor ! ExpressionEvaluationFailedResponse(this, new RuntimeException(f.toList.mkString(", ")))
    }

    WorkflowExecutionDiff(Map(this -> ExecutionStatus.Running)).validNel
  }
}

object ExpressionKey {
  private [execution] case class ExpressionEvaluationSucceededResponse(expressionKey: ExpressionKey, values: Map[OutputPort, WomValue])
  private [execution] case class ExpressionEvaluationFailedResponse(expressionKey: ExpressionKey, reason: Throwable)
}
