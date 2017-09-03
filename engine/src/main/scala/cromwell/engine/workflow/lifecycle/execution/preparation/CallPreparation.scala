package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CromwellGraphNode._
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.expression.IoFunctionSet

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object CallPreparation {
  sealed trait CallPreparationActorCommands
  case class Start(outputStore: OutputStore) extends CallPreparationActorCommands

  trait CallPreparationActorResponse

  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends CallPreparationActorResponse

  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse

  def resolveAndEvaluateInputs(callKey: CallKey,
                               workflowDescriptor: EngineWorkflowDescriptor,
                               expressionLanguageFunctions: IoFunctionSet,
                               outputStore: OutputStore): Try[Map[String, WdlValue]] = {
    import cats.instances.list._
    import cats.syntax.traverse._
    val call = callKey.scope
    val womOutputStore = outputStore.toWomOutputStore

    val inputMappingsFromPreviousCalls: Map[String, WdlValue] = call.inputPorts.map { inputPort =>
      val outputPort = inputPort.upstream
      // TODO WOM: scatters ?
      outputPort.unqualifiedName -> womOutputStore.get(outputPort).getOrElse(throw new IllegalStateException(s"Can't find output value for $outputPort"))
    } toMap
    
    // Note to self: maybe what we want is all portBasedInputs + expressionBasedInput
    // Expression based we have already, we could add portBasedInputs by finding the corresponding value of their upstream output port
    // We still have a Map[String, WdlValue] though because expressionBasedInput only gives us Strings
    // Note to self 2: looks like all inputs is actually call.callable.inputs
    // Should we map those to either expressionBasedInputs or portBased, and if no match error ? (which I think couldn't happen since wom validated that already?)
    // Should the final type be Map[InputDefinition, WdlValue] ? Will need to coerce to the InputDefinition declared type too.
    // We might lose the ability to get an FQN for the inputs then. Does it matter ?
    
    val allInputMappings = inputMappingsFromPreviousCalls ++ workflowDescriptor.backendDescriptor.knownValues
    val evaluated: Map[String, ErrorOr[WdlValue]] = call.expressionBasedInputs map {
      case (name, expr) => name -> expr.expression.evaluateValue(allInputMappings, expressionLanguageFunctions)
    }

    evaluated.toList.traverse[ErrorOr, (String, WdlValue)]({case (a, b) => b.map(c => (a,c))}).map(_.toMap) match {
      case Valid(res) => Success(res)
      case Invalid(f) => Failure(new Exception(f.toList.mkString(", ")))
    }
  }
}
