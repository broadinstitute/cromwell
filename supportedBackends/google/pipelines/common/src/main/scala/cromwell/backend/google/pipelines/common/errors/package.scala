package cromwell.backend.google.pipelines.common

package object errors {

  private def quotaMessages = List(
    "A resource limit has delayed the operation",
    "usage too high",
    "no available zones",
    "resource_exhausted",
    "quota too low",
    "waiting for quota"
  )

  def isQuotaMessage(msg: String): Boolean =
    quotaMessages.exists(msg.contains)

}
