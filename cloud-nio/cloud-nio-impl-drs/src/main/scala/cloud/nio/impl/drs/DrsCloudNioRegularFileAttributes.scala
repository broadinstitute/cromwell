package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime

import cloud.nio.spi.CloudNioRegularFileAttributes
import com.typesafe.config.Config

import scala.concurrent.duration._


class DrsCloudNioRegularFileAttributes(config: Config, drsPath: String) extends CloudNioRegularFileAttributes{

  private val drsPathResolver: DrsPathResolver = DrsPathResolver(config)

  private def runtimeException(missingKey: String) = {
    throw new RuntimeException(s"Failed to resolve DRS path $drsPath. The response from Martha doesn't contain the key '$missingKey'.")
  }

  override def fileHash: Option[String] = {
    val checksumsArrayOption = drsPathResolver.resolveDrsThroughMartha(drsPath).dos.data_object.checksums

    checksumsArrayOption match {
      case Some(checksumsArray) => checksumsArray.collectFirst{ case c if c.`type`.equalsIgnoreCase("md5") => c.checksum }
      case None => None
    }
  }

  override def lastModifiedTime(): FileTime = {
    val lastModifiedInStringOption = drsPathResolver.resolveDrsThroughMartha(drsPath).dos.data_object.updated

    lastModifiedInStringOption match {
      case Some(lastModifiedInString) => {
        val lastModifiedInDuration = Duration.apply(lastModifiedInString).toMillis
        FileTime.fromMillis(lastModifiedInDuration)
      }
      case None => runtimeException("updated")
    }
  }

  override def size(): Long = {
    val sizeOption = drsPathResolver.resolveDrsThroughMartha(drsPath).dos.data_object.size

    sizeOption match {
      case Some(size) => size
      case None => runtimeException("size")
    }
  }

  override def fileKey(): String = drsPath
}
