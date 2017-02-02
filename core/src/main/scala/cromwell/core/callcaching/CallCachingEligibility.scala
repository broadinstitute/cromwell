package cromwell.core.callcaching

sealed trait CallCachingEligibility
case object CallCachingEligible extends CallCachingEligibility
case class CallCachingIneligible(reason: String) extends CallCachingEligibility
