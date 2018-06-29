package cromwell.filesystems.gcs.bucket

/**
  * Represents knowledge of the RP value of a bucket.
  */
sealed trait RequesterPaysValue {
  /**
    * True if the project should be set, false otherwise
    */
  def withProject: Boolean

  /**
    * True for all but Disabled. This is useful when a request fails to know if it should be retried with a project or not.
    */
  def enabled: Boolean
}
object RequesterPaysValue{

  /**
    * The bucket is known to have RP on or off.
    * We can use this to confidently set (or not) a billing project in the request.
    */
  case class Known(value: Boolean) extends RequesterPaysValue {
    override def withProject = value
    override def enabled = true
  }

  /**
    * We don't know whether the bucket has RP on or not. In this case we try WITHOUT the project, and retry with it after
    */
  case object Unknown extends RequesterPaysValue {
    override def withProject = false
    override def enabled = true
  }

  /**
    * Requester pays is disabled, which means we should NEVER set the project
    */
  case object Disabled extends RequesterPaysValue {
    override def withProject = false
    override def enabled = false
  }
}
