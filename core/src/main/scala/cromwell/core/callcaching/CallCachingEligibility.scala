package cromwell.core.callcaching

sealed trait CallCachingEligibility
case object CallCachingEligible extends CallCachingEligibility
sealed trait CallCachingIneligible extends CallCachingEligibility {
  def message: String
}
  
case class FloatingDockerTagWithHash(hash: String) extends CallCachingIneligible {
  override val message = s"""You are using a floating docker tag in this task. Cromwell does not consider tasks with floating tags to be eligible for call caching.
        |If you want this task to be eligible for call caching in the future, use a docker runtime attribute with a digest instead.
        |This is the exact docker image that was used for this job: $hash
        |You can replace the docker runtime attribute in your task with the above value to make this task eligible for call caching.""".stripMargin
}

case object FloatingDockerTagWithoutHash extends CallCachingIneligible {
  override val message = s"""You are using a floating docker tag in this task. Cromwell does not consider tasks with floating tags to be eligible for call caching.
         |If you want this task to be eligible for call caching in the future, use a docker runtime attribute with a digest instead.
         |Cromwell attempted to retrieve the current hash for this docker image but failed.
         |This is not necessarily a cause for concern as Cromwell is currently only able to retrieve hashes for Dockerhub and GCR images.
         |The job will be dispatched to the appropriate backend that will attempt to run it.""".stripMargin
}
