package cromwell.core.path

import java.io.Writer

object JavaWriterImplicits {
  implicit class FlushingAndClosingWriter(writer: Writer) {
    /** Convenience method to flush and close in one shot. */
    def flushAndClose() = {
      writer.flush()
      writer.close()
    }
  }
}
