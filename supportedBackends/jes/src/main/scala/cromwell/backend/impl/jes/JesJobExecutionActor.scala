package cromwell.backend.impl.jes

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardSyncExecutionActor, StandardSyncExecutionActorParams}
import cromwell.core.Dispatcher.BackendDispatcher

/** A default implementation of the sync params. */
case class JesSyncExecutionActorParams
(
  override val jobDescriptor: BackendJobDescriptor,
  jesConfiguration: JesConfiguration,
  jesBackendInitializationData: JesBackendInitializationData,
  override val serviceRegistryActor: ActorRef,
  jesBackendSingletonActorOption: Option[ActorRef]
) extends StandardSyncExecutionActorParams {
  override val jobIdKey: String = JesJobExecutionActor.JesOperationIdKey
  override val asyncJobExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] = classOf[Nothing]
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
  override val backendInitializationDataOption: Option[BackendInitializationData] = Option(jesBackendInitializationData)
}

object JesJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            jesWorkflowInfo: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef,
            jesBackendSingletonActor: Option[ActorRef]): Props = {
    val params = JesSyncExecutionActorParams(
      jobDescriptor,
      jesWorkflowInfo,
      initializationData,
      serviceRegistryActor,
      jesBackendSingletonActor)
    Props(new JesJobExecutionActor(params)).withDispatcher(BackendDispatcher)
  }

  val JesOperationIdKey = "__jes_operation_id"
}

case class JesJobExecutionActor(jesParams: JesSyncExecutionActorParams)
  extends StandardSyncExecutionActor(jesParams) {

  override def createAsyncRefName(): String = "JesAsyncBackendJobExecutionActor"

  override def createAsyncProps(): Props = jabjeaProps

  private[jes] def jabjeaProps = {
    Props(
      new JesAsyncBackendJobExecutionActor(
        JesAsyncExecutionActorParams(
          jesParams.jobDescriptor,
          jesParams.jesConfiguration,
          jesParams.jesBackendInitializationData,
          jesParams.serviceRegistryActor,
          jesParams.jesBackendSingletonActorOption,
          completionPromise)
      )
    )
  }
}
