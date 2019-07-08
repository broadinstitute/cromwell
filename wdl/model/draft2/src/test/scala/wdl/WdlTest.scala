package wdl

import better.files.File
import better.files.File.currentWorkingDirectory
import org.scalatest.{Matchers, WordSpecLike}
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model._
import wom.ResolvedImportRecord

trait WdlTest extends Matchers with WordSpecLike {
  def resolver(root: File)(relPath: String): Draft2ResolvedImportBundle =
    Draft2ResolvedImportBundle((root / relPath).contentAsString, ResolvedImportRecord((root / relPath).pathAsString))
  def loadWdl(path: String) = loadWdlFile(currentWorkingDirectory/"wom"/"src"/"test"/"resources"/path)
  def loadWdlFile(wdlFile: File) =
    WdlNamespaceWithWorkflow.load(wdlFile.contentAsString, Seq(resolver(wdlFile / "..") _)).get
  def getTask(ns: WdlNamespace, name: String): WdlTask = ns.tasks.find(_.unqualifiedName == name).get
  def getCall(ns: WdlNamespaceWithWorkflow, name: String): WdlTaskCall = ns.workflow.taskCalls.find(_.unqualifiedName == name) getOrElse {
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
