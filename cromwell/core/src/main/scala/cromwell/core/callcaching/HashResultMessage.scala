package cromwell.core.callcaching


object HashKey {
  def apply(keyComponents: String*) = new HashKey(true, keyComponents.toList)
  def apply(checkForHitOrMiss: Boolean, keyComponents: String*) = new HashKey(checkForHitOrMiss, keyComponents.toList)
}

case class HashKey(checkForHitOrMiss: Boolean, keyComponents: List[String]) {
  val key = keyComponents.mkString(": ")
}
case class HashValue(value: String)
case class HashResult(hashKey: HashKey, hashValue: HashValue)

sealed trait HashResultMessage
trait SuccessfulHashResultMessage extends HashResultMessage {
  def hashes: Set[HashResult]
}
case class MultiHashingFailedMessage(keys: Set[HashKey], reason: Throwable) extends HashResultMessage
case class HashingFailedMessage(file: String, reason: Throwable) extends HashResultMessage
case object HashingServiceUnvailable extends HashResultMessage
