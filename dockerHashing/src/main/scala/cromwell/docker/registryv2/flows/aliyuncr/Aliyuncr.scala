package cromwell.docker.registryv2.flows.aliyuncr

import scala.util.matching.Regex

object AliyunCr {
  
  val supportAliunCrRegion = List("cn-qingdao", "cn-beijing", "cn-zhangjiakou", "cn-huhehaote",
    "cn-hangzhou", "cn-shanghai", "cn-shenzhen", "cn-hongkong",
    "ap-northeast-1", "ap-southeast-1", "ap-southeast-2", "ap-southeast-3",
    "ap-southeast-5", "ap-south-1", "us-east-1", "us-west-1", "me-east-1",
    "eu-central-1").mkString("|")
  val validAliyunCrHosts: Regex = ("""registry.""" + s"""(?:($supportAliunCrRegion))""" + """.aliyuncs.com""").r

  def isValidAliyunCrHost(host: Option[String]): Boolean =
    host match {
      case Some(h) => { 
        h match {
          case validAliyunCrHosts(_) => true
          case _ => false
        }
      }
      case _ => false
    }
}
