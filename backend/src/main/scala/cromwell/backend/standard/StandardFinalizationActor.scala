package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.core.CallOutputs
import wdl4s.TaskCall

trait StandardFinalizationActorParams {
  def workflowDescriptor: BackendWorkflowDescriptor

  def calls: Set[TaskCall]
  
  def serviceRegistryActor: ActorRef

  def jobExecutionMap: JobExecutionMap

  def workflowOutputs: CallOutputs

  def initializationDataOption: Option[BackendInitializationData]

  def configurationDescriptor: BackendConfigurationDescriptor
}

case class DefaultStandardFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  serviceRegistryActor: ActorRef,
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData],
  configurationDescriptor: BackendConfigurationDescriptor
) extends StandardFinalizationActorParams

trait StandardFinalizationActor extends BackendWorkflowFinalizationActor {

  def standardParams: StandardFinalizationActorParams

  override lazy val workflowDescriptor: BackendWorkflowDescriptor = standardParams.workflowDescriptor
  override lazy val calls: Set[TaskCall] = standardParams.calls
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor
  lazy val initializationDataOption: Option[BackendInitializationData] = standardParams.initializationDataOption
  lazy val jobExecutionMap: JobExecutionMap = standardParams.jobExecutionMap
  lazy val workflowOutputs: CallOutputs = standardParams.workflowOutputs
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor
}
