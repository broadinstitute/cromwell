package cromwell.backend.impl.spark

import java.nio.file.Path

import better.files._
import com.typesafe.scalalogging.StrictLogging
import cromwell.core.path.PathImplicits._
import cromwell.core.path.{TailedWriter, UntailedWriter}

import scala.sys.process._
import scala.util.{Failure, Success, Try}

object SparkCommands {
  val SparkSubmit = "spark-submit"
  val Separator = "--"
  val ExecutorCores = "total-executor-cores"
  val ExecutorMemory = "executor-memory"
  val AppMainClass = "class"
  val Master = "master"
  val DeployMode = "deploy-mode"
  val SparkAppWithArgs = "spark_app_with_args"
  val MemoryUnit = "G" //Accepts only GB
}

class SparkCommands extends StrictLogging {
  import SparkCommands._
  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  def writeScript(instantiatedCommand: String, filePath: Path, containerRoot: Path) = {

    val scriptBody =
      s"""

#!/bin/sh
cd $containerRoot
$instantiatedCommand
echo $$? > rc

       """.trim + "\n"
    File(filePath).write(scriptBody)
  }

  def sparkSubmitCommand(attributes: Map[String, Any]): String = {
    val sparkHome = Try(sys.env("SPARK_HOME")) match {
      case Success(s) => Option(s)
      case Failure(ex) =>
        logger.warn(s"Spark home does not exist picking up default command")
        None
    }

    val sparkSubmit = sparkHome
      .map(p => s"$p/bin/$SparkSubmit")
      .getOrElse(SparkSubmit)

    attributes.foreach {
      case (key, _) => attributes.get(key) match {
        case Some(v) => s" $Separator$key $v "
        case None =>
          logger.warn (s" key: $key doesn't exist spark should not pick it")
      }
    }

    val mapWithoutSparkAppKey = attributes - SparkAppWithArgs

    val partialSparkCmd = mapWithoutSparkAppKey.foldLeft(s"$sparkSubmit "){
      (acc, pair) => s"$acc $Separator${pair._1} ${mapWithoutSparkAppKey(pair._1)} "
    }

    val cmdWithSparkApp = s"$partialSparkCmd ${attributes(SparkAppWithArgs)}"

    logger.debug(s"spark command is : $cmdWithSparkApp")
    cmdWithSparkApp
  }
}

trait SparkProcess {

  def processStdout: String = stdout.toString.trim

  def processStderr: String = stderr.toString.trim

  def commandList(command: String): Seq[String] = List("/bin/bash", command)

  def untailedWriter(path: Path): UntailedWriter = path.untailed

  def tailedWriter(limit: Int, path: Path): TailedWriter = path.tailed(limit)

  def externalProcess(cmdList: Seq[String], processLogger: ProcessLogger): Process = cmdList.run(processLogger)
}
