package cromwell.engine.backend.jes

/**
 * Trait for JesConfiguration
 */
trait JesConfiguration {
  def jesConf: JesAttributes
}

object ProductionJesConfiguration {
  lazy val jesConf =  JesAttributes()
}

trait ProductionJesConfiguration extends JesConfiguration {
  override lazy val jesConf = ProductionJesConfiguration.jesConf
}