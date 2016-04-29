package cromwell.engine.backend
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleTaskAbortedException extends Exception {
  override def getMessage: String = "The task was aborted."
}
