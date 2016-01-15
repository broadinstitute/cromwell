package cromwell.engine.backend

import cromwell.engine.io.IoInterface
import cromwell.engine.{CallContext, CallEngineFunctions, WorkflowContext, WorkflowEngineFunctions}

import scala.language.postfixOps

class DefaultWorkflowEngineFunctions(interface: IoInterface, context: WorkflowContext) extends WorkflowEngineFunctions(interface, context)

class DefaultCallEngineFunctions(interface: IoInterface, context: CallContext) extends DefaultWorkflowEngineFunctions(interface ,context)
with CallEngineFunctions
