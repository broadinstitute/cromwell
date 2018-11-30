package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime

import cloud.nio.spi.CloudNioRegularFileAttributes
import com.typesafe.config.Config

import scala.concurrent.duration._


class DrsCloudNioRegularFileAttributes(config: Config, drsPath: String) extends CloudNioRegularFileAttributes{

  private val drsPathResolver: DrsPathResolver = DrsPathResolver(config)

  override def fileHash: Option[String] = {
    val checksumsArray = drsPathResolver.resolveDrsThroughMartha(drsPath).dos.data_object.checksums

    checksumsArray.collectFirst{ case c if c.`type`.equalsIgnoreCase("md5") => c.checksum }
  }

  override def lastModifiedTime(): FileTime = {
    val lastModifiedInString = drsPathResolver.resolveDrsThroughMartha(drsPath).dos.data_object.updated
    val lastModifiedInDuration = Duration.apply(lastModifiedInString).toMillis

    FileTime.fromMillis(lastModifiedInDuration)
  }

  override def size(): Long = drsPathResolver.resolveDrsThroughMartha(drsPath).dos.data_object.size

  override def fileKey(): String = drsPath
}
