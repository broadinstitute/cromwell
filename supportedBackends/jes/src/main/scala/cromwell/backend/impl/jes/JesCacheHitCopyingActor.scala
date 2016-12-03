package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import cromwell.backend.callcaching.CacheHitDuplicating
import cromwell.backend.{BackendCacheHitCopyingActor, BackendJobDescriptor}
import cromwell.core.path.PathCopier
import cromwell.core.logging.JobLogging

case class JesCacheHitCopyingActor(override val jobDescriptor: BackendJobDescriptor,
                                   jesConfiguration: JesConfiguration,
                                   initializationData: JesBackendInitializationData,
                                   serviceRegistryActor: ActorRef)
  extends BackendCacheHitCopyingActor with CacheHitDuplicating with JesJobCachingActorHelper with JobLogging {
  override protected def duplicate(source: Path, destination: Path) = PathCopier.copy(source, destination).get

  override protected def destinationCallRootPath = jesCallPaths.callExecutionRoot

  override protected def destinationJobDetritusPaths = jesCallPaths.detritusPaths
  
  override val workflowDescriptor = jobDescriptor.workflowDescriptor
}

object JesCacheHitCopyingActor {

  def props(jobDescriptor: BackendJobDescriptor,
            jesConfiguration: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesCacheHitCopyingActor(jobDescriptor, jesConfiguration, initializationData, serviceRegistryActor))
  }
}
