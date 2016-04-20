package cromwell.engine.backend.dummy

import akka.actor.Props
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionSucceededResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.{BackendJobExecutionAbortSucceededResponse, JobAbortResponse}
import cromwell.backend._
import cromwell.core.CallOutput
import wdl4s.types._
import wdl4s.values._
import wdl4s.TaskOutput

import scala.concurrent.Future

object DummyBackendJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = Props(DummyBackendJobExecutionActor(jobDescriptor, configurationDescriptor))
}

case class DummyBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor, override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {
  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    val outputs = (jobDescriptor.call.task.outputs map taskOutputToJobOutput).toMap
    Future.successful(BackendJobExecutionSucceededResponse(jobDescriptor.key, outputs))
  }

  private def taskOutputToJobOutput(taskOutput: TaskOutput) =
    taskOutput.name -> CallOutput(sampleValue(taskOutput.wdlType), None)

  private def sampleValue(wdlType: WdlType): WdlValue = wdlType match {
    case WdlIntegerType => WdlInteger(3)
    case WdlFloatType => WdlFloat(55.55)
    case WdlStringType => WdlString("The rain in Spain falls mainly in the plain")
    case WdlBooleanType => WdlBoolean(true)
    case WdlFileType => WdlFile("/root/of/all/evil")
    case WdlArrayType(memberType) => WdlArray(WdlArrayType(memberType), List(sampleValue(memberType)))
    case WdlObjectType => WdlObject(Map("a" -> WdlString("1"), "b" -> WdlString("2")))
    case WdlMapType(keyType, valueType) => WdlMap(WdlMapType(keyType, valueType), Map(sampleValue(keyType) -> sampleValue(valueType)))
  }

  /**
    * Restart or resume a previously-started job.
    */
  override def recover = execute

  /**
    * Abort a running job.
    */
  override def abortJob: Future[JobAbortResponse] = Future.successful(BackendJobExecutionAbortSucceededResponse(jobDescriptor.key))

}
