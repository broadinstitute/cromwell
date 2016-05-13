package cromwell.core

import java.sql.Timestamp
import java.util.Date


trait KnowsWhatTimeItIs {
  def now = new Timestamp(new Date().getTime)
}
