package cromwell.filesystems.gcs

import java.nio.file.OpenOption

object ContentTypeOption {
  sealed trait ContentType
  case object PlainText extends ContentType with OpenOption {
    override def toString = "plain/text"
  }
  case object Json extends ContentType with OpenOption {
    override def toString = "application/json"
  }
}


