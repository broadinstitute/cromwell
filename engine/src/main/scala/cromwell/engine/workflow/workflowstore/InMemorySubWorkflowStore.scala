package cromwell.engine.workflow.workflowstore

import cromwell.database.sql.tables.SubWorkflowStoreEntry
import cromwell.subworkflowstore.SubWorkflowStore

import scala.concurrent.{ExecutionContext, Future}

class InMemorySubWorkflowStore(workflowStore: InMemoryWorkflowStore) extends SubWorkflowStore {

  var subWorkflowStore: Set[SubWorkflowStoreEntry] = Set.empty

  override def addSubWorkflowStoreEntry(rootWorkflowExecutionUuid: String,
                                        parentWorkflowExecutionUuid: String,
                                        callFullyQualifiedName: String,
                                        jobIndex: Int,
                                        jobAttempt: Int,
                                        subWorkflowExecutionUuid: String
  )(implicit ec: ExecutionContext): Future[Unit] =
    if (workflowStore.workflowStore.exists { case (wf, _) => wf.id.toString == rootWorkflowExecutionUuid }) {
      subWorkflowStore = subWorkflowStore +
        SubWorkflowStoreEntry(
          rootWorkflowId = Option(0),
          parentWorkflowExecutionUuid = parentWorkflowExecutionUuid,
          callFullyQualifiedName = callFullyQualifiedName,
          callIndex = jobIndex,
          callAttempt = jobAttempt,
          subWorkflowExecutionUuid = subWorkflowExecutionUuid
        )
      Future.successful(())
    } else Future.failed(new Throwable(s"No such root workflow: $rootWorkflowExecutionUuid"))

  override def querySubWorkflowStore(parentWorkflowExecutionUuid: String,
                                     callFqn: String,
                                     jobIndex: Int,
                                     jobAttempt: Int
  )(implicit ec: ExecutionContext): Future[Option[SubWorkflowStoreEntry]] =
    Future.successful(
      subWorkflowStore.find(k =>
        k.parentWorkflowExecutionUuid == parentWorkflowExecutionUuid &&
          k.callFullyQualifiedName == callFqn &&
          k.callIndex == jobIndex &&
          k.callAttempt == jobAttempt
      )
    )

  override def removeSubWorkflowStoreEntries(
    parentWorkflowExecutionUuid: String
  )(implicit ec: ExecutionContext): Future[Int] = {
    val toRemove = subWorkflowStore.filter(k => k.parentWorkflowExecutionUuid == parentWorkflowExecutionUuid)
    subWorkflowStore = subWorkflowStore -- toRemove
    Future.successful(toRemove.size)
  }
}
