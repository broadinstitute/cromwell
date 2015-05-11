package cromwell.engine.backend.local

import java.io.{BufferedWriter, File, FileWriter}

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, TaskOutput}
import cromwell.engine.backend.Backend
import cromwell.engine.SymbolStore

import scala.language.postfixOps
import scala.util.Try

class LocalBackend extends Backend {

  override def executeCommand(commandLine: String, call: Call, taskOutputs: Set[TaskOutput], symbolStore: SymbolStore): Map[String, Try[WdlValue]] = {
    import sys.process._

    // TODO make these proper temp files, but for now it's preferable
    // TODO to leave these in a deterministic place.
    val stdoutFile = new File("stdout.txt")
    val stderrFile = new File("stderr.txt")
    val stdout = new BufferedWriter(new FileWriter(stdoutFile))
    val stderr = new BufferedWriter(new FileWriter(stderrFile))

    commandLine ! ProcessLogger(stdout write, stderr write)

    Seq(stdout, stderr).foreach { writer =>
      writer.flush()
      writer.close()
    }

    taskOutputs.map { taskOutput =>
      taskOutput.name -> taskOutput.expression.evaluate(
        scopedLookupFunction(call, symbolStore),
        new LocalEngineFunctions(new TaskExecutionContext(stdoutFile, stderrFile, new File("."))))
    }.toMap
  }
}
