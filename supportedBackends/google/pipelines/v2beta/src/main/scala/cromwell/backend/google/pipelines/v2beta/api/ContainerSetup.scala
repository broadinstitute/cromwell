package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk

trait ContainerSetup {
  def containerSetupActions(mounts: List[Mount]): List[Action] = {
    val containerRoot = PipelinesApiWorkingDisk.MountPoint.pathAsString

    // As opposed to V1, the container root does not have a 777 umask, which can cause issues for docker running as non root
    // Run a first action to create the root and give it the right permissions
    val containerRootSetup = ActionBuilder
      .cloudSdkShellAction(s"mkdir -p $containerRoot && chmod -R a+rwx $containerRoot")(
        mounts = mounts,
        labels = Map(Key.Tag -> Value.ContainerSetup)
      )

    ActionBuilder.annotateTimestampedActions("container setup", Value.ContainerSetup)(List(containerRootSetup))
  }
}
