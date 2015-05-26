package cromwell.singleuser

import java.io.File

import cromwell.Main
import cromwell.binding.WdlBinding
import cromwell.server.WorkflowManagerSystem

import scala.util.{Failure, Success, Try}

// FIXME: Once we have a definition of what 'run' means in Main.scala beyond some pretty printing, it goes here
object SingleWorkflowRunner extends WorkflowManagerSystem {

  // FIXME: This function is only here temporarily, it'll never really exist in this form (or this name)
  def runIt(args: Array[String]): Unit = {
    if (args.isEmpty) Main.usageAndExit()
    else {
      Try(WdlBinding.process(new File(args(1)))) match {
        case Success(b) => Main.printProcessedBinding(b)
        case Failure(e) =>
          println(e)
          System.exit(-1)
      }
    }
  }
}
