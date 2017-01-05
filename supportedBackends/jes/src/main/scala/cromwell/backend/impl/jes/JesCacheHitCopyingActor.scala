package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.backend.callcaching.CacheHitDuplicating
import cromwell.backend.{BackendCacheHitCopyingActor, BackendConfigurationDescriptor, BackendJobDescriptor, BackendWorkflowDescriptor}
import cromwell.core.path.PathCopier
import cromwell.core.logging.JobLogging

case class JesCacheHitCopyingActor(override val jobDescriptor: BackendJobDescriptor,
                                   jesConfiguration: JesConfiguration,
                                   initializationData: JesBackendInitializationData,
                                   serviceRegistryActor: ActorRef)
  extends BackendCacheHitCopyingActor with CacheHitDuplicating with JesJobCachingActorHelper with JobLogging {
  override protected def duplicate(source: Path, destination: Path): Unit = PathCopier.copy(source, destination).get

  override protected lazy val destinationCallRootPath: Path = jesCallPaths.callExecutionRoot

  override protected lazy val destinationJobDetritusPaths: Map[String, Path] = jesCallPaths.detritusPaths
  
  override val workflowDescriptor: BackendWorkflowDescriptor = jobDescriptor.workflowDescriptor

  override lazy val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

object JesCacheHitCopyingActor {

  def props(jobDescriptor: BackendJobDescriptor,
            jesConfiguration: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesCacheHitCopyingActor(jobDescriptor, jesConfiguration, initializationData, serviceRegistryActor)).withDispatcher(BackendDispatcher)
  }
}
