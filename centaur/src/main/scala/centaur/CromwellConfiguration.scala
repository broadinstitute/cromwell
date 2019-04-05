package centaur

trait CromwellProcess {
  def logFile: String
  def displayString: String
  def start(): Unit
  def stop(): Unit
  def isAlive: Boolean
  def cromwellConfiguration: CromwellConfiguration
}

trait CromwellConfiguration {
  def createProcess: CromwellProcess
  def logFile: String
}
