package centaur

import better.files.File
import centaur.api.CentaurCromwellClient

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.sys.process._

/**
  * This is not thread-safe, and hence only works if cromwell is started / stopped from only 1 thread at a time.
  */
object CromwellManager {
  val ManagedCromwellPort = 8008
  private val timeout = 30 seconds
  private val interval = 1 second
  private var cromwellProcess: Option[Process] = None
  private var processLogger: FileProcessLogger = _
  private var _ready: Boolean = false
  private var _isManaged: Boolean = false
  
  /**
    * Returns true if Cromwell is ready to be queried, false otherwise
    * In Unmanaged mode, this is irrelevant so always return true.
    * In managed mode return the value of _ready
    */
  def isReady: Boolean = !_isManaged || _ready
  
  // Check that we have a cromwellProcess, that this process is alive, and that cromwell is ready to accept requests 
  private def isAlive = cromwellProcess.exists(_.isAlive()) && CentaurCromwellClient.isAlive

  def startCromwell(cromwell: CromwellConfiguration): Unit = {
    _isManaged = true
    
    if (!isAlive) {
      val cromwellCommand = s"java -Dconfig.file=${cromwell.conf} -Dwebservice.port=$ManagedCromwellPort -jar ${cromwell.jar} server"
      
      // Redirect stdout and stderr to the log file
      val logFile: File = File(cromwell.logFile)
      processLogger = ProcessLogger(logFile.toJava)
        
      // Start the cromwell process
      println("Starting Cromwell...")
      val process = cromwellCommand.run(processLogger)
      cromwellProcess = Option(process)

      var waitedFor = Duration.Zero

      while (!isAlive && waitedFor < timeout) {
        println("Waiting for Cromwell...")
        Thread.sleep(interval.toMillis)
        waitedFor = waitedFor + interval
      }
      
      _ready = true
      if (isAlive) println("Cromwell is running")
      else {
        println("Timeout waiting for cromwell server - failing test run")
        println(logFile.contentAsString)
        stopCromwell()
        System.exit(1)
      }
    }
  }

  def stopCromwell() = {
    _ready = false
    println("Stopping Cromwell...")
    cromwellProcess foreach { process =>
      processLogger.flush()
      processLogger.close()
      process.destroy()
      process.exitValue()
    }
    
    cromwellProcess = None
  }
}
