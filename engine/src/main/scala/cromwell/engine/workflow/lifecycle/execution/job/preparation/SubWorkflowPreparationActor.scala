package cromwell.engine.workflow.lifecycle.execution.job.preparation

import akka.actor.{Actor, Props}
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import common.exception.MessageAggregation
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobBreadCrumb
import cromwell.core.Dispatcher._
import cromwell.core.WorkflowId
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation._
import cromwell.engine.workflow.lifecycle.execution.job.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import cromwell.engine.workflow.lifecycle.execution.keys.SubWorkflowKey
import cromwell.engine.{EngineIoFunctions, EngineWorkflowDescriptor}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{GraphInputNode, GraphNode, OptionalGraphInputNode, OptionalGraphInputNodeWithDefault}
import wom.values.{WomEvaluatedCallInputs, WomValue}

class SubWorkflowPreparationActor(workflowDescriptor: EngineWorkflowDescriptor,
                                  expressionLanguageFunctions: EngineIoFunctions,
                                  callKey: SubWorkflowKey,
                                  subWorkflowId: WorkflowId
) extends Actor
    with WorkflowLogging {

  override lazy val workflowIdForLogging = workflowDescriptor.possiblyNotRootWorkflowId
  override lazy val rootWorkflowIdForLogging = workflowDescriptor.rootWorkflowId

  private def prepareExecutionActor(inputEvaluation: WomEvaluatedCallInputs): ErrorOr[CallPreparationActorResponse] = {
    val oldBackendDescriptor = workflowDescriptor.backendDescriptor

    evaluateStartingKnownValues(inputEvaluation, callKey.node.callable.graph.inputNodes) map { startingValues =>
      val newBackendDescriptor = oldBackendDescriptor.copy(
        id = subWorkflowId,
        callable = callKey.node.callable,
        knownValues = startingValues,
        breadCrumbs = oldBackendDescriptor.breadCrumbs :+ BackendJobBreadCrumb(workflowDescriptor.callable,
                                                                               workflowDescriptor.id,
                                                                               callKey
        )
      )
      val engineDescriptor =
        workflowDescriptor.copy(backendDescriptor = newBackendDescriptor, parentWorkflow = Option(workflowDescriptor))
      SubWorkflowPreparationSucceeded(engineDescriptor, inputEvaluation)
    }
  }

  /**
    * Work out a set of  "workflow inputs" to pass to this subworkflow as though it were a top-level workflow receiving inputs
    */
  private def evaluateStartingKnownValues(inputEvaluation: WomEvaluatedCallInputs,
                                          workflowInputs: Set[GraphInputNode]
  ): ErrorOr[Map[OutputPort, WomValue]] = {

    // Find the values in the provided inputs that match up with subworkflow inputs. Silently drop the rest on the floor.
    val providedInputs: Map[OutputPort, WomValue] = inputEvaluation.toList.flatMap { case (name, value) =>
      workflowInputs.collectFirst {
        case i if i.singleOutputPort.identifier.localName.value == name.localName.value => i.singleOutputPort -> value
      }
    }.toMap

    // We don't need to supply these values, but we might have. Used once, in the validation step below.
    def optionalsAndDefaults: Set[GraphNode] = workflowInputs collect {
      case optionalInputNode: OptionalGraphInputNode => optionalInputNode
      case optionalWithDefault: OptionalGraphInputNodeWithDefault => optionalWithDefault
    }

    // Make sure the subworkflow will be getting a value for every required input:
    NonEmptyList.fromList(
      (workflowInputs.toSet[GraphNode] diff providedInputs.keySet.map(_.graphNode) diff optionalsAndDefaults).toList
    ) match {
      case None => Valid(providedInputs)
      case Some(missingNodeNel) =>
        Invalid(
          missingNodeNel map (n => s"Couldn't find starting value for subworkflow input: ${n.identifier.localName}")
        )
    }
  }

  override def receive = {
    case Start(valueStore) =>
      val evaluatedInputs = resolveAndEvaluateInputs(callKey, expressionLanguageFunctions, valueStore)
      import common.validation.ErrorOr._
      evaluatedInputs.flatMap(prepareExecutionActor) match {
        case Valid(response) => context.parent ! response
        case Invalid(f) =>
          context.parent ! CallPreparationFailed(
            callKey,
            new MessageAggregation {
              override def exceptionContext: String = "Failed to evaluate inputs for sub workflow"
              override def errorMessages: Iterable[String] = f.toList
            }
          )
      }
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }
}

object SubWorkflowPreparationActor {
  case class SubWorkflowPreparationSucceeded(workflowDescriptor: EngineWorkflowDescriptor,
                                             inputs: WomEvaluatedCallInputs
  ) extends CallPreparationActorResponse

  def props(workflowDescriptor: EngineWorkflowDescriptor,
            expressionLanguageFunctions: EngineIoFunctions,
            key: SubWorkflowKey,
            subWorkflowId: WorkflowId
  ) =
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new SubWorkflowPreparationActor(workflowDescriptor, expressionLanguageFunctions, key, subWorkflowId))
      .withDispatcher(EngineDispatcher)
}
