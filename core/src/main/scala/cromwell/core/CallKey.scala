package cromwell.core

import wdl4s.Call

trait CallKey extends JobKey {
  def scope: Call
}
