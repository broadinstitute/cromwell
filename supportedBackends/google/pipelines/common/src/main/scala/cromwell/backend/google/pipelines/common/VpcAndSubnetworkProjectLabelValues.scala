package cromwell.backend.google.pipelines.common

final case class VpcAndSubnetworkProjectLabelValues(projectId: Option[String], region:Option[String], vpcName: String, subnetNameOpt: Option[String])
