package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.Mount
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.ToParameter.ops._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.Value

trait Localization {
  def localizeActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration) = {
    val jobInputLocalization = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toActions(mounts).toList)
    ActionBuilder.annotateTimestampedActions("localization", Value.Localization)(jobInputLocalization)
  }
}
