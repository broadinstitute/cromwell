package cromwell.backend.google.pipelines.common

final case class VpcAndSubnetworkProjectLabelValues(sharedProjectId: Option[String], sharedRegion:Option[String], vpcName: String, subnetNameOpt: Option[String])
