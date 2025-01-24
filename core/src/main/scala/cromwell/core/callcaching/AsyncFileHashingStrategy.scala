package cromwell.core.callcaching

sealed trait AsyncFileHashingStrategy

object AsyncFileHashingStrategy {
  case object Crc32c extends AsyncFileHashingStrategy
  case object Md5 extends AsyncFileHashingStrategy
  case object Md5ThenIdentity extends AsyncFileHashingStrategy
  case object ETag extends AsyncFileHashingStrategy

  // TODO validate fs type here?
  def apply(s: String): Option[AsyncFileHashingStrategy] = s.toLowerCase() match {
    case "md5" => Some(Md5)
    case "crc32c" => Some(Crc32c)
    case "md5+identity" => Some(Md5ThenIdentity)
    case "etag" => Some(ETag)
    case _ => None
  }
}
