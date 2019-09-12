package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait MemoryRetryCheckAction {

  def checkForMemoryRetryActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    createPipelineParameters.runtimeAttributes.retryWithDoubleMemory match {
      case Some(keys) => List(ActionBuilder.checkForMemoryRetryAction(keys, mounts))
      case None => List.empty[Action]
    }
  }
}
