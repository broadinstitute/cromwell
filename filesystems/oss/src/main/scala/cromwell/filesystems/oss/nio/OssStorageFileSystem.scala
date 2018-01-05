package cromwell.filesystems.oss.nio


import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.util.Objects
import java.{lang, util}

import com.aliyun.oss.OSSClient

import scala.collection.JavaConverters._


object OssStorageFileSystem {
  val SEPARATOR: String = "/"
  val URI_SCHEMA: String = "oss"
  val OSS_VIEW = "oss"
  val BASIC_VIEW = "basic"

  def apply(provider: OssStorageFileSystemProvider, bucket: String, config: OssStorageConfiguration): OssStorageFileSystem = {
    val res = new OssStorageFileSystem(bucket, config)
    res.internalProvider = provider
    res
  }
}

object OssStorageConfiguration {
  val ENDPOINT_KEY = "endpoint"
  val ACCESS_ID_KEY = "access-id"
  val ACCESS_KEY_KEY = "access-key"
  val SECURITY_TOKEN_KEY = "security-token"

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

    new OssStorageConfiguration(endpoint, accessId, accessKey, securityToken)
  }

  def getClient(map: Map[String, String]): OSSClient = {
    parseMap(map).newOssClient()
  }

  def getClient(endpoint: String,
                accessId: String,
                accessKey: String,
                stsToken: Option[String]): OSSClient = {
    OssStorageConfiguration(endpoint, accessId, accessKey, stsToken).newOssClient()
  }

}

final case class OssStorageConfiguration(endpoint: String,
                                   accessId: String,
                                   accessKey: String,
                                   stsToken: Option[String] = None
                                  ) {
  import OssStorageConfiguration._

  def toMap: Map[String, String] = {
    val ret = Map(ENDPOINT_KEY -> endpoint, ACCESS_ID_KEY -> accessId, ACCESS_KEY_KEY -> accessKey)
    val token = stsToken map {token => SECURITY_TOKEN_KEY -> token}
    ret ++ token
  }

  def newOssClient() = {
    stsToken match {
      case Some(token: String) =>
        new OSSClient(endpoint, accessId, accessKey, token)
      case None =>
        new OSSClient(endpoint, accessId, accessKey)
    }
  }
}

case class OssStorageFileSystem(bucket: String, config: OssStorageConfiguration) extends FileSystem {

  var internalProvider: OssStorageFileSystemProvider = OssStorageFileSystemProvider(config)

  override def provider: OssStorageFileSystemProvider = internalProvider

  override def getPath(first: String, more: String*): OssStoragePath = OssStoragePath.getPath(this, first, more: _*)

  override def close(): Unit = {
    // do nothing currently.
  }

  override def isOpen: Boolean = {
    true
  }

  override def isReadOnly: Boolean = {
    false
  }

  override def getSeparator: String = {
    OssStorageFileSystem.SEPARATOR
  }

  override def getRootDirectories: lang.Iterable[Path] = {
    Set[Path](OssStoragePath.getPath(this, UnixPath.ROOT_PATH)).asJava
  }

  override def getFileStores: lang.Iterable[FileStore] = {
    Set.empty[FileStore].asJava
  }

  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = {
    FileSystems.getDefault.getPathMatcher(syntaxAndPattern)
  }

  override def getUserPrincipalLookupService: UserPrincipalLookupService = {
    throw new UnsupportedOperationException()
  }

  override def newWatchService(): WatchService = {
    throw new UnsupportedOperationException()
  }

  override def supportedFileAttributeViews(): util.Set[String] = {
    Set(OssStorageFileSystem.OSS_VIEW, OssStorageFileSystem.BASIC_VIEW).asJava
  }

  override def equals(obj: scala.Any): Boolean = {
    this == obj ||
      obj.isInstanceOf[OssStorageFileSystem] &&
        obj.asInstanceOf[OssStorageFileSystem].config.equals(config) &&
        obj.asInstanceOf[OssStorageFileSystem].bucket.equals(bucket)
  }

  override def hashCode(): Int = Objects.hash(bucket)

  override def toString: String = OssStorageFileSystem.URI_SCHEMA + "://" + bucket
}

