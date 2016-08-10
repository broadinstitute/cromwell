package cromwell.core.callcaching

case class HashKey(key: String, checkForHitOrMiss: Boolean = true)
case class HashValue(value: String)
case class HashResult(hashKey: HashKey, hashValue: HashValue)

sealed trait HashResultMessage
trait SuccessfulHashResultMessage extends HashResultMessage {
  def hashes: Set[HashResult]
}
case class HashingFailedMessage(key: HashKey, reason: Throwable) extends HashResultMessage