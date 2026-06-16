package cromwell.backend.impl.tes

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import net.ceedubs.ficus.Ficus._
import org.apache.commons.codec.digest.DigestUtils

import java.nio.file.{Files, Paths}
import scala.util.{Failure, Try, Using}

object TesBackendFileHashingActor {
  def props(standardParams: StandardFileHashingActorParams): Props =
    Props(new TesBackendFileHashingActor(standardParams))

  /**
    * Returns true for an absolute, non-cloud path that belongs to the shared filesystem
    * inside a TES worker container.
    *
    * When `localRoot` is provided (from `filesystems.local.local-root` config), only
    * paths under that specific mount point are considered local. This avoids accidentally
    * treating any absolute path (e.g. `/etc/hostname`) as a shared-filesystem file.
    *
    * When `localRoot` is None (default), falls back to the previous generic behaviour:
    * any path starting with "/" that has no URI scheme is considered local.
    * This preserves backward compatibility with deployments that do not set local-root.
    */
  def isLocalPath(value: String, localRoot: Option[String] = None): Boolean = {
    val isAbsoluteLocal = value.startsWith("/") && !value.startsWith("//") && !value.contains("://")
    localRoot match {
      case Some(root) => isAbsoluteLocal && value.startsWith(root)
      case None       => isAbsoluteLocal
    }
  }
}

/**
  * Custom file-hashing actor for the TES backend.
  *
  * For local-style paths (absolute paths without a URI scheme, e.g. /mnt/efs/…):
  *   - If `filesystems.local.caching.check-sibling-md5 = true` and a sibling
  *     `<file>.md5` exists, the content of that file is returned as the hash.
  *   - Otherwise the file is hashed with MD5.
  *
  * For cloud / HTTP / DRS / … paths the method returns `None` so the base-class
  * async IO hashing takes over (S3 ETags, GCS CRC32, …).
  */
class TesBackendFileHashingActor(standardParams: StandardFileHashingActorParams)
    extends StandardFileHashingActor(standardParams) {

  // Use the default IoCommandBuilder – cloud paths are handled by the super-class.
  override val ioCommandBuilder = cromwell.core.io.DefaultIoCommandBuilder

  /** Read check-sibling-md5 from the backend configuration. */
  lazy val checkSiblingMd5: Boolean =
    configurationDescriptor.backendConfig
      .as[Option[Config]]("filesystems.local.caching")
      .flatMap(_.as[Option[Boolean]]("check-sibling-md5"))
      .getOrElse(false)

  /**
    * Mount point of the shared filesystem inside TES worker containers.
    * Reads `filesystems.local.local-root` (preferred) with backward-compat fallback
    * to legacy key `filesystems.local.efs`.
    * When set, only paths under this root are hashed locally; other absolute paths are
    * treated as cloud paths and handed to the base-class async IO hashing.
    * PR note: replaces the implicit /mnt/efs convention.
    */
  lazy val localRoot: Option[String] =
    configurationDescriptor.backendConfig
      .as[Option[String]]("filesystems.local.local-root")
      .orElse(configurationDescriptor.backendConfig.as[Option[String]]("filesystems.local.efs"))

  /**
    * Called before async IO hashing.  Returning `Some(result)` short-circuits
    * the async path; returning `None` falls through to it.
    */
  override def customHashStrategy(fileRequest: SingleFileHashRequest): Option[Try[String]] = {
    val rawValue = fileRequest.file.valueString
    if (TesBackendFileHashingActor.isLocalPath(rawValue, localRoot)) {
      Some(hashLocalFile(rawValue))
    } else {
      // Cloud / HTTP / DRS paths → let the base class handle them asynchronously.
      None
    }
  }

  // ---------------------------------------------------------------------------
  // Private helpers
  // ---------------------------------------------------------------------------

  private def hashLocalFile(pathStr: String): Try[String] = {
    val nioPath = Paths.get(pathStr)

    if (!Files.exists(nioPath))
      return Failure(new java.io.FileNotFoundException(s"Local file not found for hashing: $pathStr"))

    if (checkSiblingMd5) {
      val md5SiblingPath = Paths.get(pathStr + ".md5")
      if (Files.exists(md5SiblingPath)) {
        log.debug(s"[TES-hash] Using sibling .md5 file for: $pathStr")
        Try(new String(Files.readAllBytes(md5SiblingPath)).trim)
      } else {
        log.debug(s"[TES-hash] No sibling .md5 found, computing MD5 for: $pathStr")
        Using(Files.newInputStream(nioPath))(DigestUtils.md5Hex)
      }
    } else {
      log.debug(s"[TES-hash] Computing MD5 for: $pathStr")
      Using(Files.newInputStream(nioPath))(DigestUtils.md5Hex)
    }
  }
}
