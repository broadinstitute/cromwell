package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait EgressCheckAction {
  def egressCheckAction(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    Nil
  }
}