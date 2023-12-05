package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.BackendInitializationData
import cromwell.backend.io.{JobPaths, WorkflowPaths}

import scala.concurrent.ExecutionContext

class StandardInitializationData(
  val workflowPaths: WorkflowPaths,
  val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  val standardExpressionFunctionsClass: Class[_ <: StandardExpressionFunctions]
) extends BackendInitializationData {

  /* TODO: This could (should?) be monadic instead of reflection. */
  private lazy val standardExpressionFunctionsConstructor =
    standardExpressionFunctionsClass.getConstructor(classOf[StandardExpressionFunctionsParams])

  def expressionFunctions(jobPaths: JobPaths,
                          ioActorProxy: ActorRef,
                          ec: ExecutionContext
  ): StandardExpressionFunctions = {
    val pathBuilders = jobPaths.workflowPaths.pathBuilders
    val standardParams = DefaultStandardExpressionFunctionsParams(pathBuilders, jobPaths.callContext, ioActorProxy, ec)
    standardExpressionFunctionsConstructor.newInstance(standardParams)
  }
}
