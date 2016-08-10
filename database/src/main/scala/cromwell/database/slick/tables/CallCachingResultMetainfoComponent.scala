package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingResultMetaInfoEntry

trait CallCachingResultMetaInfoComponent {

  this: DriverComponent =>

  import driver.api._

  class CallCachingResultMetaInfoEntries(tag: Tag) extends Table[CallCachingResultMetaInfoEntry](tag, "CALL_CACHING_RESULT_METAINFO") {
    def callCachingResultMetaInfoId = column[Int]("CALL_CACHING_RESULT_METAINFO_ID", O.PrimaryKey, O.AutoInc)
    def workflowUuid = column[String]("WORKFLOW_UUID")
    def callFqn = column[String]("CALL_FQN")
    def returnCode = column[Option[Int]]("RETURN_CODE")
    def scatterIndex = column[Int]("JOB_SCATTER_INDEX")
    def allowResultReuse = column[Boolean]("ALLOW_RESULT_REUSE")

    override def * = (workflowUuid, callFqn, scatterIndex, returnCode, allowResultReuse, callCachingResultMetaInfoId.?) <>
      (CallCachingResultMetaInfoEntry.tupled, CallCachingResultMetaInfoEntry.unapply)

    def ccmUniquenessConstraint = index("UK_CALL_CACHING_RESULT_METAINFO", (workflowUuid, callFqn, scatterIndex), unique = true)
  }

  protected val callCachingResultMetaInfo = TableQuery[CallCachingResultMetaInfoEntries]

  val callCachingResultMetaInfoAutoInc = callCachingResultMetaInfo returning callCachingResultMetaInfo.
    map(_.callCachingResultMetaInfoId) into ((a, id) => a.copy(callCachingResultMetaInfoEntryId = Option(id)))

  /**
    * Useful for finding the call caching result meta info for a given ID
    */
  val metaInfoById = Compiled(
    (metaInfoId: Rep[Int]) => for {
      metaInfo <- callCachingResultMetaInfo
      if metaInfo.callCachingResultMetaInfoId === metaInfoId
    } yield metaInfo)
}
