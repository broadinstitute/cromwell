package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import cromwell.backend.google.batch.io.GcpBatchWorkingDisk

trait ContainerSetup {
  import RunnableLabels._

  def containerSetupRunnables(volumes: List[Volume]): List[Runnable] = {
    val containerRoot = GcpBatchWorkingDisk.MountPoint.pathAsString

    // As opposed to V1, the container root does not have a 777 umask, which can cause issues for docker running as non root
    // Run a first action to create the root and give it the right permissions
    val containerRootSetup = RunnableBuilder
      .cloudSdkShellRunnable(s"mkdir -p $containerRoot && chmod -R a+rwx $containerRoot")(
        volumes = volumes,
        labels = Map(Key.Tag -> Value.ContainerSetup),
        flags = List.empty
      )

    RunnableBuilder
      .annotateTimestampedRunnable("container setup", Value.ContainerSetup)(volumes, List(containerRootSetup))
      .map(_.build())
  }
}
