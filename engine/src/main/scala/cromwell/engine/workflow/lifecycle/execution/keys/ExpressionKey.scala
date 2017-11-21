package cromwell.engine.workflow.lifecycle.execution.keys

import akka.actor.ActorRef
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import cromwell.core.ExecutionIndex._
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.keys.ExpressionKey.{ExpressionEvaluationFailedResponse, ExpressionEvaluationSucceededResponse}
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActorData, WorkflowExecutionDiff}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExpressionNode
import wom.values.WomValue

/**
  * Key for expression nodes.
  */
private [execution] case class ExpressionKey(node: ExpressionNode, index: ExecutionIndex) extends JobKey {
  override val attempt = 1
  override lazy val tag = s"Expression-${node.localName}:${index.fromIndex}:$attempt"

  val womType = node.womType
  val singleOutputPort: OutputPort = node.singleExpressionOutputPort

  lazy val upstreamPorts: Map[String, OutputPort] = node.inputMapping map {
    case (key, input) => key -> input.upstream
  }

  def processRunnable(data: WorkflowExecutionActorData, workflowExecutionActor: ActorRef): ErrorOr[WorkflowExecutionDiff] = {
    upstreamPorts.traverseValues(data.valueStore.resolve(index)) map { lookup =>
      // Send a message to self in case we decide to change evaluate to return asynchronously, if we don't we could
      // directly add the value to the value store in the execution diff
      node.evaluateAndCoerce(lookup, data.expressionLanguageFunctions) match {
        case Right(result) => workflowExecutionActor ! ExpressionEvaluationSucceededResponse(this, result)
        case Left(f) => workflowExecutionActor ! ExpressionEvaluationFailedResponse(this, new RuntimeException(f.toList.mkString(", ")))
      }
    } valueOr { f =>
      workflowExecutionActor ! ExpressionEvaluationFailedResponse(this, new RuntimeException(s"Unable to start $this: " + f.toList.mkString(", ")))
    }

    WorkflowExecutionDiff(Map(this -> ExecutionStatus.Running)).validNel
  }
}

object ExpressionKey {
  private [execution] case class ExpressionEvaluationSucceededResponse(declarationKey: ExpressionKey, value: WomValue)
  private [execution] case class ExpressionEvaluationFailedResponse(declarationKey: ExpressionKey, reason: Throwable)
}
