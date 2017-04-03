package cromwell.core.callcaching

sealed trait MaybeCallCachingEligible {
  def dockerHash: Option[String]
}

sealed trait CallCachingEligible extends MaybeCallCachingEligible
sealed trait CallCachingIneligible extends MaybeCallCachingEligible

case object NoDocker extends CallCachingEligible {
  override def dockerHash: Option[String] = None
}
case class DockerWithHash(dockerAttribute: String) extends CallCachingEligible {
  override def dockerHash: Option[String] = Option(dockerAttribute)
}

case class FloatingDockerTagWithoutHash(dockerTag: String) extends CallCachingIneligible {
  override def dockerHash: Option[String] = None
}
