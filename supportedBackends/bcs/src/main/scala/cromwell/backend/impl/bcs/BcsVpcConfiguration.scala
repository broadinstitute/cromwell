package cromwell.backend.impl.bcs


import scala.util.{Failure, Success, Try}


final case class BcsVpcConfiguration(cidrBlock: Option[String] = None,
                               vpcId: Option[String] = None)


object BcsVpcConfiguration {
  val cidrBlockPattern = """(([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\.){3}([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])(\/([0-9]|[1-2][0-9]|3[0-2]){1,2})""".r
  val vpcIdPattern = """(vpc-[^\s]+)""".r

  def parse(s: String): Try[BcsVpcConfiguration] = {
    val cidrBlock = cidrBlockPattern findFirstIn s
    val vpcId = vpcIdPattern findFirstIn s

    if (cidrBlock.isEmpty && vpcId.isEmpty) {
      Failure(new IllegalArgumentException("vpc configuration must be a string like '192.168.0.0/16 vpc-xxxx' "))
    } else {
      Success(BcsVpcConfiguration(cidrBlock, vpcId))
    }
  }
}