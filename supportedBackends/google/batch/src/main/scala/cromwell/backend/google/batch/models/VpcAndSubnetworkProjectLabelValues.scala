package cromwell.backend.google.batch.models

import cromwell.backend.google.batch.models.VpcAndSubnetworkProjectLabelValues._

final case class VpcAndSubnetworkProjectLabelValues(vpcName: String, subnetNameOpt: Option[String]) {

  /**
   * Returns a qualified network name replacing the string `\${projectId}` in the network name if found.
   */
  def networkName(projectId: String): String = {
    val networkNameTemplate =
      if (vpcName.contains("/")) {
        vpcName
      } else {
        s"projects/$ProjectIdToken/global/networks/$vpcName/"
      }

    networkNameTemplate.replace(ProjectIdToken, projectId)
  }

  /**
   * Replaces the string `\${projectId}` in the subnet name if found.
   */
  def subnetNameOption(projectId: String): Option[String] =
    subnetNameOpt map { _.replace(ProjectIdToken, projectId) }
}

object VpcAndSubnetworkProjectLabelValues {
  private val ProjectIdToken = s"$${projectId}"
}
