package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder.Labels.{Key, Value}

import scala.collection.JavaConverters._

trait ContainerSetup {
  def containerSetupActions(mounts: List[Mount]): List[Action] = {
    val containerRoot = PipelinesApiWorkingDisk.MountPoint.pathAsString

    // As opposed to V1, the container root does not have a 777 umask, which can cause issues for docker running as non root
    // Run a first action to create the root and give it the right permissions
    val containerRootSetup = ActionBuilder
      .cloudSdkAction
      .setCommands(
        List("/bin/bash", "-c", s"mkdir -p $containerRoot && chmod -R a+rwx $containerRoot").asJava
      )
      .setMounts(mounts.asJava)
      .setLabels(Map(Key.Tag -> Value.ContainerSetup).asJava)

    ActionBuilder.annotateTimestampedActions("container setup", Value.ContainerSetup)(List(containerRootSetup))
  }
}
