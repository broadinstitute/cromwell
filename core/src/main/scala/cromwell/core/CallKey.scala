package cromwell.core

import wdl4s.wdl.WdlCall

trait CallKey extends JobKey {
  def scope: WdlCall
}
