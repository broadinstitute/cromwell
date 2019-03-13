package centaur

import better.files.File
import centaur.api.CentaurCromwellClient
import com.typesafe.scalalogging.StrictLogging

import scala.concurrent.duration._
import scala.language.postfixOps

/**
  * This is not thread-safe, and hence only works if cromwell is started / stopped from only 1 thread at a time.
  */
object CromwellManager extends StrictLogging {
  val ManagedCromwellPort = 8008
  val timeout = 120 seconds
  private val interval = 5 second
  private val timeoutExitStatus = 66
  private var cromwellProcess: Option[CromwellProcess] = None
  private var _ready: Boolean = false
  private var _isManaged: Boolean = false
  
  /**
    * Returns true if Cromwell is ready to be queried, false otherwise
    * In Unmanaged mode, this is irrelevant so always return true.
    * In managed mode return the value of _ready
    */
  def isReady: Boolean = !_isManaged || _ready
  
  // Check that we have a cromwellProcess, that this process is alive, and that cromwell is ready to accept requests 
  private def isAlive(checkType: String): Boolean = {
    val processAlive = cromwellProcess.exists(_.isAlive)
    logger.info(s"Cromwell process alive $checkType = $processAlive")
    if (processAlive) {
      val serverAlive = CentaurCromwellClient.isAlive
      logger.info(s"Cromwell server alive $checkType = $serverAlive")
      serverAlive
    } else {
      false
    }
  }

  def startCromwell(cromwellConfiguration: CromwellConfiguration): Unit = {
    _isManaged = true

    if (!isAlive("at start")) {
      cromwellProcess = Option(cromwellConfiguration.createProcess)
      cromwellProcess foreach { _.start() }

      var waitedFor = Duration.Zero
      var seenAlive = false

      def wasOrIsAlive(): Boolean = {
        if (!seenAlive && isAlive("while waiting")) {
          seenAlive = true
        }
        seenAlive
      }

      while (!wasOrIsAlive() && waitedFor < timeout) {
        logger.info("Waiting for Cromwell...")
        Thread.sleep(interval.toMillis)
        waitedFor = waitedFor + interval
      }

      _ready = true
      if (wasOrIsAlive()) logger.info("Cromwell is running")
      else {
        logger.error("Timeout waiting for cromwell server - failing test run")
        logger.error(File(cromwellConfiguration.logFile).contentAsString)
        stopCromwell("Timed out waiting for server")
        System.exit(timeoutExitStatus)
      }
    }
  }

  def stopCromwell(reason: String) = {
    _ready = false
    logger.info(s"Stopping Cromwell... ($reason)")
    try {
      cromwellProcess foreach { _.stop() }
    } catch {
      case e: Exception => 
        logger.error("Caught exception while stopping Cromwell")
        e.printStackTrace()
    }
    
    cromwellProcess = None
  }
}
