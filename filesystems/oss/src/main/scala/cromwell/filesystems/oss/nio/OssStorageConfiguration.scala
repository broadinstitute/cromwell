package cromwell.filesystems.oss.nio

import com.aliyun.oss.{ClientConfiguration, OSSClient}
import com.aliyun.oss.common.auth.DefaultCredentialProvider

object OssStorageConfiguration {
  val ENDPOINT_KEY = "endpoint"
  val ACCESS_ID_KEY = "access-id"
  val ACCESS_KEY_KEY = "access-key"
  val SECURITY_TOKEN_KEY = "security-token"
  val REFRESH_INTERVAL = "refresh-interval"

  import scala.collection.immutable.Map
  def parseMap(map: Map[String, Any]): OssStorageConfiguration = {
    val endpoint = map.get(ENDPOINT_KEY) match {
      case Some(endpoint: String) if !endpoint.isEmpty => endpoint
      case _ => throw new IllegalArgumentException(s"endpoint is mandatory and must be an unempty string")
    }
    val accessId = map.get(ACCESS_ID_KEY) match {
      case Some(id: String) if !id.isEmpty => id
      case _ => throw new IllegalArgumentException(s"access-id is mandatory and must be an unempty string")
    }

    val accessKey = map.get(ACCESS_KEY_KEY) match {
      case Some(key: String) if !key.isEmpty => key
      case _ => throw new IllegalArgumentException(s"access-key is mandatory and must be an unempty string")
    }

    val securityToken = map.get(SECURITY_TOKEN_KEY) match {
      case Some(token: String) if !token.isEmpty => Some(token)
      case _ => None
    }

    new DefaultOssStorageConfiguration(endpoint, accessId, accessKey, securityToken)
  }

  def getClient(map: Map[String, String]): OSSClient = {
    parseMap(map).newOssClient()
  }

  def getClient(endpoint: String,
                accessId: String,
                accessKey: String,
                stsToken: Option[String]): OSSClient = {
    DefaultOssStorageConfiguration(endpoint, accessId, accessKey, stsToken).newOssClient()
  }

}

trait OssStorageConfiguration {

  import OssStorageConfiguration._

  def endpoint: String
  def accessId: String
  def accessKey: String
  def stsToken: Option[String]

  def toMap: Map[String, String] = {
    val ret = Map(
      ENDPOINT_KEY -> endpoint,
      ACCESS_ID_KEY -> accessId,
      ACCESS_KEY_KEY -> accessKey)

    val token = stsToken map {token => SECURITY_TOKEN_KEY -> token}
    ret ++ token
  }

  def newOssClient() = {
    stsToken match {
      case Some(token: String) =>
        val cred = new DefaultCredentialProvider(accessId, accessKey, token)

        new OSSClient(endpoint, cred, new ClientConfiguration)
      case None =>
        val cred = new DefaultCredentialProvider(accessId, accessKey)
        new OSSClient(endpoint, cred, new ClientConfiguration)
    }
  }
}

final case class DefaultOssStorageConfiguration(val endpoint: String, val accessId: String, val accessKey: String, val stsToken: Option[String]= None) extends OssStorageConfiguration

