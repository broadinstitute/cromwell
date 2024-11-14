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
        s"projects/$ProjectIdToken/global/networks/$vpcName"
      }

    networkNameTemplate.replace(ProjectIdToken, projectId)
  }

  /**
   * Replaces the string `\${projectId}` in the subnet name if found.  Replace wildcard character used in terra configuration for subnetworks with appropriate region
   */
  def subnetNameOption(projectId: String, region: String): Option[String] = {
    subnetNameOpt map { _.replace(ProjectIdToken, projectId) }
    subnetNameOpt map { _.replace(regionWildcard, "regions/"+ region) }
  }
}

object VpcAndSubnetworkProjectLabelValues {
  private val ProjectIdToken = s"$${projectId}"
  private val regionWildcard = s"regions/*"
}
