package cromwell.core

import java.sql.Timestamp
import java.util.Date


trait KnowsWhatTimeItIs {
  def currentTime = new Timestamp(new Date().getTime)
}
