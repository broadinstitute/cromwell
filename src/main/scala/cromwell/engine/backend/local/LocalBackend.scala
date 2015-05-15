package cromwell.engine.backend.local

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, TaskOutput}
import cromwell.engine.SymbolStore
import cromwell.engine.backend.Backend
import cromwell.util.FileUtil
import cromwell.util.FileUtil._

import scala.language.postfixOps
import scala.sys.process._
import scala.util.Try

class LocalBackend extends Backend {

  override def executeCommand(commandLine: String, call: Call, taskOutputs: Set[TaskOutput], symbolStore: SymbolStore): Map[String, Try[WdlValue]] = {

    val (stdoutFile, stdoutWriter) = FileUtil.tempFileAndWriter("stdout")
    val (stderrFile, stderrWriter) = FileUtil.tempFileAndWriter("stderr")
    val (commandFile, commandWriter) = FileUtil.tempFileAndWriter("command")

    commandWriter.write(commandLine)
    commandWriter.flushAndClose()

    s"/bin/bash $commandFile" ! ProcessLogger(stdoutWriter write, stderrWriter write)

    Seq(stdoutWriter, stderrWriter).foreach {
      _.flushAndClose()
    }

    taskOutputs.map { taskOutput =>
      taskOutput.name -> taskOutput.expression.evaluate(
        scopedLookupFunction(call, symbolStore),
        new LocalEngineFunctions(TaskExecutionContext(stdoutFile, stderrFile)))
    }.toMap
  }
}
