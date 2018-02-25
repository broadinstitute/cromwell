package centaur.cwl

import centaur.CentaurConfig

object CentaurCwlRunnerConfig {
  lazy val conf = CentaurConfig.conf.getConfig("cwl-runner")
}
