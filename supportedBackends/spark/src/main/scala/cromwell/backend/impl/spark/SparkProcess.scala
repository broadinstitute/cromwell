package cromwell.backend.impl.spark

import java.nio.file.Path
import com.typesafe.scalalogging.StrictLogging
import cromwell.core.{TailedWriter, UntailedWriter}
import cromwell.core.PathFactory.EnhancedPath
import scala.sys.process._
import scala.language.postfixOps

class SparkProcess extends StrictLogging {
  private val stdout = new StringBuilder
  private val stderr = new StringBuilder

  def processLogger: ProcessLogger = ProcessLogger(stdout append _, stderr append _)
  def processStdout: String = stdout.toString().trim
  def processStderr: String = stderr.toString().trim
  def commandList(command: String): Seq[String] = Seq("/bin/bash",command)
  def untailedWriter(path: Path): UntailedWriter = path.untailed
  def tailedWriter(limit: Int, path: Path): TailedWriter = path.tailed(limit)
  def externalProcess(cmdList: Seq[String], processLogger: ProcessLogger = processLogger): Process = cmdList.run(processLogger)
}
