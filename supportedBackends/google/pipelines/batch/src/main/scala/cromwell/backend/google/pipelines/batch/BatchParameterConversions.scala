package cromwell.backend.google.pipelines.batch

import cromwell.backend.google.pipelines.batch.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import com.google.cloud.batch.v1.Volume
import simulacrum.typeclass

@typeclass trait ToParameter[A <: BatchParameter] {
  def toRunnables(p: A, mounts: List[Volume])(implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Volume]
}

trait BatchParameterConversions {



}
