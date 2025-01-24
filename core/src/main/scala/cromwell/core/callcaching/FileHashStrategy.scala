package cromwell.core.callcaching

// File hashing strategies used by IoHashCommand, primarily when obtaining file hashes
// for call caching purposes.
sealed trait FileHashStrategy

object FileHashStrategy {
  case object Crc32c extends FileHashStrategy
  case object Md5 extends FileHashStrategy
  case object Md5ThenIdentity extends FileHashStrategy
  case object ETag extends FileHashStrategy

  // TODO validate fs type here?
  def apply(s: String): Option[FileHashStrategy] = s.toLowerCase() match {
    case "md5" => Some(Md5)
    case "crc32c" => Some(Crc32c)
    case "md5+identity" => Some(Md5ThenIdentity)
    case "etag" => Some(ETag)
    case _ => None
  }
}
