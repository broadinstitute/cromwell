package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait MemoryRetryCheckAction {

  def checkForMemoryRetryActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    createPipelineParameters.retryWithMoreMemoryKeys match {
      case Some(keys) => List(ActionBuilder.checkForMemoryRetryAction(keys, mounts))
      case None => List.empty[Action]
    }
  }
}
