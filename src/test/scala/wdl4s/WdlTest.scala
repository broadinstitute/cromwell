package wdl4s

import better.files.File
import org.scalatest.{Matchers, WordSpecLike}

trait WdlTest extends Matchers with WordSpecLike {
  def resolver(root: File)(relPath: String): String = (root / relPath).contentAsString
  def loadWdlFile(wdlFile: File) = WdlNamespaceWithWorkflow.load(wdlFile.path, resolver(wdlFile / "..") _)
  def getTask(ns: WdlNamespace, name: String): Task = ns.tasks.find(_.unqualifiedName == name).get
  def getCall(ns: WdlNamespaceWithWorkflow, name: String): TaskCall = ns.workflow.taskCalls.find(_.unqualifiedName == name) getOrElse {
    fail(s"Expecting call with name '$name'")
  }
  def getScatter(ns: WdlNamespaceWithWorkflow, index: Int): Scatter = {
    val scatterFqn = ns.workflow.unqualifiedName + ".$scatter_" + index
    val resolution = ns.resolve(scatterFqn).collect({ case s: Scatter => s})
    resolution getOrElse {
      fail(s"Expecting a scatter block with FQN $scatterFqn")
    }
  }
  def getIf(ns: WdlNamespaceWithWorkflow, index: Int): If = {
    val ifFqn = ns.workflow.unqualifiedName + ".$if_" + index
    val resolution = ns.resolve(ifFqn).collect({ case i: If => i })
    resolution getOrElse {
      fail(s"Expecting a scatter block with FQN $ifFqn")
    }
  }
}
