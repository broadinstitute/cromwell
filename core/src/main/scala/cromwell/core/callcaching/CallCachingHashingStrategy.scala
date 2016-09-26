package cromwell.core.callcaching

sealed trait CallCachingHashingStrategy

case class ExistingSiblingHash(extension: String) extends CallCachingHashingStrategy
case object PathHash extends CallCachingHashingStrategy
case object ComputeHash extends CallCachingHashingStrategy

