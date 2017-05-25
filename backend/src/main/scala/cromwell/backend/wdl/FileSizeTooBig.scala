package cromwell.backend.wdl

case class FileSizeTooBig(override val getMessage: String) extends Exception

