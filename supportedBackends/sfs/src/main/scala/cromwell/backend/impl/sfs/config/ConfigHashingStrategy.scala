package cromwell.backend.impl.sfs.config

import java.io.{FileNotFoundException, InputStream}

import akka.event.LoggingAdapter
import com.typesafe.config.Config
import cromwell.backend.standard.StandardInitializationData
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.path.{Path, PathFactory}
import cromwell.util.TryWithResource._
import net.ceedubs.ficus.Ficus._
import net.jpountz.xxhash.XXHashFactory
import org.apache.commons.codec.digest.DigestUtils
import org.slf4j.{Logger, LoggerFactory}

import scala.util.{Failure, Try}

object ConfigHashingStrategy {
  val logger: Logger = LoggerFactory.getLogger(getClass)
  val defaultStrategy = HashFileMd5Strategy(checkSiblingMd5 = false)

  def apply(hashingConfig: Config): ConfigHashingStrategy = {
      val checkSiblingMd5 = hashingConfig.as[Option[Boolean]]("check-sibling-md5").getOrElse(false)

      // Fingerprint strategy by default checks the first 10 MiB (10485760 bytes) for performance reasons.
      // 100 MB will take to much time on network file systems. 1 MB might not be unique enough.
      // The value is user configurable.
      lazy val fingerprintSize = hashingConfig.as[Option[Long]]("fingerprint-size").getOrElse(10L * 1024 * 1024)

      hashingConfig.as[Option[String]]("hashing-strategy").getOrElse("file") match {
        case "path" => HashPathStrategy(checkSiblingMd5)
        case "file" => HashFileMd5Strategy(checkSiblingMd5)
        case "md5" => HashFileMd5Strategy(checkSiblingMd5)
        case "path+modtime" => HashPathModTimeStrategy(checkSiblingMd5)
        case "xxh64" => HashFileXxH64Strategy(checkSiblingMd5)
        case "fingerprint" => FingerprintStrategy(checkSiblingMd5, fingerprintSize)
        case what =>
          logger.warn(s"Unrecognized hashing strategy $what.")
          HashPathStrategy(checkSiblingMd5)
      }
  }
}

abstract class ConfigHashingStrategy {
  def checkSiblingMd5: Boolean
  protected def hash(file: Path): Try[String]
  protected def description: String

  protected lazy val checkSiblingMessage: String =
    if (checkSiblingMd5) "Check first for sibling md5 and if not found " else ""

  def getHash(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    def usingStandardInitData(initData: StandardInitializationData) = {
      val pathBuilders = initData.workflowPaths.pathBuilders
      val file = PathFactory.buildPath(request.file.valueString, pathBuilders).followSymbolicLinks
      if (!file.exists) Failure(new FileNotFoundException(s"Cannot hash file $file because it can't be found")) else {
        if (checkSiblingMd5) {
          precomputedMd5(file) match {
            case Some(md5) => Try(md5.contentAsString.trim)
            case None => hash(file)
          }
        } else hash(file)
      }
    }

    request.initializationData match {
      case Some(initData: StandardInitializationData) => usingStandardInitData(initData)
      case _ => Failure(new IllegalArgumentException("Need SharedFileSystemBackendInitializationData to calculate hash."))
    }
  }

  private def precomputedMd5(file: Path): Option[Path] = {
    val md5 = file.sibling(s"${file.name}.md5")
    if (md5.exists) Option(md5) else None
  }

  override def toString: String = {
    s"Call caching hashing strategy: $checkSiblingMessage$description."
  }
}

final case class HashPathStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: Path): Try[String] = {
    Try(DigestUtils.md5Hex(file.toAbsolutePath.pathAsString))
  }

  override val description = "hash file path"
}

final case class HashPathModTimeStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: Path): Try[String] = {
    // Add the last modified date here to make sure these are the files we are looking for.
    Try(DigestUtils.md5Hex(file.toAbsolutePath.pathAsString + file.lastModifiedTime.toString))
  }

  override val description = "hash file path and last modified time"
}

final case class HashFileMd5Strategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override protected def hash(file: Path): Try[String] = {
    tryWithResource(() => file.newInputStream) { DigestUtils.md5Hex }
  }

  override val description = "hash file content with md5"
}

final case class HashFileXxH64Strategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override protected def hash(file: Path): Try[String] = {
    tryWithResource(() => file.newInputStream) {HashFileXxH64StrategyMethods.xxh64sum(_)}
  }
  override val description = "hash file content with xxh64"
}

final case class FingerprintStrategy(checkSiblingMd5: Boolean, fingerprintSize: Long) extends ConfigHashingStrategy {
  override protected def hash(file: Path): Try[String] = {
    Try {
      // Calculate the xxh64 hash of last modified time and filesize. These are NOT added, as it will lead to loss of
      // information. Instead their hexstrings are concatenated and then hashed.
      HashFileXxH64StrategyMethods.xxh64sumString(file.lastModifiedTime.toEpochMilli.toHexString +
      file.size.toHexString) +
      HashFileXxH64StrategyMethods.xxh64sum(file.newInputStream, maxSize = fingerprintSize)
      }
    }
  override val description = "fingerprint the file with last modified time, size and a xxh64 hash of the first part of the file"
}

object HashFileXxH64StrategyMethods {
  // For more information about the choice of buffer size: https://github.com/rhpvorderman/hashtest/
  private lazy val defaultBufferSize: Int = 128 * 1024
  private lazy val xxhashFactory: XXHashFactory = XXHashFactory.fastestInstance()

  /**
    * Returns the xxh64sum of an input stream. The input stream is read in a buffered way.
    * @param inputStream an input Stream
    * @param bufferSize the size in bytes for the buffer.
    * @param maxSize, only calculate the hash for the first maxSize bytes. Must be a multiply of bufferSize.
    * @return A hex string of the digest.
    */
  def xxh64sum(inputStream: InputStream,
               bufferSize: Int = defaultBufferSize,
               maxSize: Long = Long.MaxValue,
               seed: Long = 0L): String = {
    val hasher = xxhashFactory.newStreamingHash64(seed)
    val buffer: Array[Byte] = new Array[Byte](bufferSize)
    var byteCounter: Long = 0
    try {
      while (inputStream.available() > 0 && byteCounter < maxSize) {
        val length: Int = inputStream.read(buffer)
        hasher.update(buffer, 0, length)
        byteCounter += length
      }
    }
    finally inputStream.close()
    // Long.toHexString does not add leading zero's
    f"%%16s".format(hasher.getValue.toHexString).replace(" ", "0")
  }

  // Only instantiate the xxh64hasher once
  private lazy val xxh64hasher = xxhashFactory.hash64()

  def xxh64sumString(string: String, seed: Long = 0L): String = {
    val bytes: Array[Byte] = string.toCharArray.map(_.toByte)
    val hash = xxh64hasher.hash(bytes, 0, bytes.length, seed)
    // Long.toHexString does not add leading zero's
    f"%%16s".format(hash.toHexString).replace(" ", "0")
  }
}

