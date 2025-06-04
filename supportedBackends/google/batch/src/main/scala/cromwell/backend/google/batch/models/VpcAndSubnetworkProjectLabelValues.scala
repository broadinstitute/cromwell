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
   * If a subnet name is found, it returns the fully qualified subnet name by replacing \${projectId} and
   * substituting any wildcard characters used in Terra's subnetwork configuration with the appropriate region.
   */
  def subnetNameOption(projectId: String, region: String): Option[String] =
    subnetNameOpt map { subnetName =>
      val subnetworkNameTemplate = if (subnetName.contains("/")) {
        subnetName
      } else {
        s"projects/$ProjectIdToken/regions/*/subnetworks/$subnetName"
      }

      subnetworkNameTemplate.replace(ProjectIdToken, projectId).replace(regionWildcard, "regions/" + region)
    }
}

object VpcAndSubnetworkProjectLabelValues {
  private val ProjectIdToken = s"$${projectId}"
  private val regionWildcard = s"regions/*"
}
