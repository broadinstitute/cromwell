package cromwell.engine.backend

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, WdlStandardLibraryFunctions, WorkflowDescriptor}

import scala.util.Try

trait BackendCall {
  val backend: Backend
  val workflowDescriptor: WorkflowDescriptor
  val call: Call
  def engineFunctions: WdlStandardLibraryFunctions
  def lookupFunction: String => WdlValue
  def instantiateCommand: Try[String]
  def execute: Try[Map[String, WdlValue]]
}
