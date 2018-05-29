package cromwell.backend

case class FileSizeTooBig(override val getMessage: String) extends Exception

