package cromwell.engine.backend.local

import java.io.Writer

import cromwell.binding.types.WdlFileType
import cromwell.binding.values.{WdlFile, WdlString, WdlValue}
import cromwell.binding.{Call, TaskOutput}
import cromwell.engine.SymbolStore
import cromwell.engine.backend.Backend
import cromwell.util.FileUtil
import cromwell.util.FileUtil._

import scala.language.postfixOps
import scala.sys.process._
import scala.util.Try


object LocalBackend {
  implicit class WriteWithNewline(val writer: Writer) extends AnyVal {
    def writeWithNewline(string: String): Unit = {
      writer.write(string)
      writer.write("\n")
    }
  }
}


class LocalBackend extends Backend {

  /**
   * Executes the specified command line, using the supplied call and symbol store for expression evaluation.
   * Returns a `Map[String, Try[WdlValue]]` of output names to values.
   */
  override def executeCommand(commandLine: String, call: Call, taskOutputs: Set[TaskOutput], symbolStore: SymbolStore): Map[String, Try[WdlValue]] = {
    import LocalBackend._

    val (stdoutFile, stdoutWriter) = FileUtil.tempFileAndWriter("stdout")
    val (stderrFile, stderrWriter) = FileUtil.tempFileAndWriter("stderr")
    val (commandFile, commandWriter) = FileUtil.tempFileAndWriter("command")

    commandWriter.write(commandLine)
    commandWriter.flushAndClose()

    s"/bin/bash $commandFile" ! ProcessLogger(stdoutWriter writeWithNewline, stderrWriter writeWithNewline)

    Vector(stdoutWriter, stderrWriter).foreach {
      _.flushAndClose()
    }

    /**
     * Handle possible auto-conversion from string literal to WdlFile
     * and possible magic "stdout" or "stderr" literals.
     */
    def possiblyAutoConvertedValue(taskOutput: TaskOutput, wdlValue: WdlValue): WdlValue = {
      wdlValue match {
        case v: WdlString =>
          // Handle possible auto-conversion from string literal to WdlFile
          // and possible magic "stdout" or "stderr" literals.
          taskOutput.wdlType match {
            case WdlFileType =>
              v.value match {
                case "stdout" => WdlFile(stdoutFile)
                case "stderr" => WdlFile(stderrFile)
                case _ => WdlFile(v.value)
              }
            case _ => v
          }
        case v => v
      }
    }

    taskOutputs.map { taskOutput =>
      val rawValue = taskOutput.expression.evaluate(
        scopedLookupFunction(call, symbolStore),
        new LocalEngineFunctions(TaskExecutionContext(stdoutFile, stderrFile)))

      taskOutput.name -> rawValue.map { possiblyAutoConvertedValue(taskOutput, _) }
    }.toMap
  }
}
