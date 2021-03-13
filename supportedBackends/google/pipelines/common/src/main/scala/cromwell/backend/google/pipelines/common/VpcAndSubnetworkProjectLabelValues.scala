package cromwell.backend.google.pipelines.common

final case class VpcAndSubnetworkProjectLabelValues(project:Option[String], vpcName: String, subnetNameOpt: Option[String])
