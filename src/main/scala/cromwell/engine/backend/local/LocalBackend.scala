package cromwell.engine.backend.local

import java.io.{Writer, BufferedWriter, File, FileWriter}

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, TaskOutput}
import cromwell.engine.backend.Backend
import cromwell.engine.SymbolStore

import scala.language.postfixOps
import scala.util.Try

class LocalBackend extends Backend {

  override def executeCommand(commandLine: String, call: Call, taskOutputs: Set[TaskOutput], symbolStore: SymbolStore): Map[String, Try[WdlValue]] = {

    implicit class FlushingAndClosingWriter(writer: Writer) {
      /** Convenience method to flush and close in one shot. */
      def flushAndClose() = {
        writer.flush()
        writer.close()
      }
    }

    def buildFileAndWriter(filename: String): (File, Writer) = {
      val file = new File(filename)
      val writer = new BufferedWriter(new FileWriter(file))
      (file, writer)
    }

    import sys.process._

    // TODO make these proper temp files, but for now it's preferable
    // TODO to leave them in a deterministic place.
    val (stdoutFile, stdoutWriter) = buildFileAndWriter("stdout.txt")
    val (stderrFile, stderrWriter) = buildFileAndWriter("stderr.txt")
    val (commandFile, commandWriter) = buildFileAndWriter("command.txt")

    commandWriter.write(commandLine)
    commandWriter.flushAndClose()

    s"/bin/bash ${commandFile.getAbsolutePath}" ! ProcessLogger(stdoutWriter write, stderrWriter write)

    Seq(stdoutWriter, stderrWriter).foreach { _.flushAndClose() }

    taskOutputs.map { taskOutput =>
      taskOutput.name -> taskOutput.expression.evaluate(
        scopedLookupFunction(call, symbolStore),
        new LocalEngineFunctions(new TaskExecutionContext(stdoutFile, stderrFile, new File("."))))
    }.toMap
  }
}
