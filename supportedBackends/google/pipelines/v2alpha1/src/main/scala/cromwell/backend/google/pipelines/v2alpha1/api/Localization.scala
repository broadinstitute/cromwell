package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.Mount
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.ToParameter.ops._

import scala.collection.JavaConverters._

trait Localization {
  def localizeActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]) = {
    val containerRoot = PipelinesApiWorkingDisk.MountPoint.pathAsString

    // As opposed to V1, the container root does not have a 777 umask, which can cause issues for docker running as non root
    // Run a first action to create the root and give it the right permissions
    val containerRootSetup = ActionBuilder
      .cloudSdkAction
      .setCommands(
        List("/bin/bash", "-c", s"mkdir -p $containerRoot && chmod -R a+rwx $containerRoot").asJava
      )
      .setMounts(mounts.asJava)

    val jobInputLocalization = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toActions(mounts, createPipelineParameters.projectId).toList)
    
    List(containerRootSetup) ++ jobInputLocalization
  }
}
