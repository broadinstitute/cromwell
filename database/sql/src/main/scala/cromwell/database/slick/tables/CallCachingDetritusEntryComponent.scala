package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingJobDetritusEntry

trait CallCachingJobDetritusComponent {

  this: DriverComponent with CallCachingResultMetaInfoComponent =>

  import driver.api._

  class CallCachingJobDetritus(tag: Tag) extends Table[CallCachingJobDetritusEntry](tag, "CALL_CACHING_JOB_DETRITUS") {
    def callCachingJobDetritusId = column[Int]("CALL_CACHING_JOB_DETRITUS_ID", O.PrimaryKey, O.AutoInc)
    def jobDetritusKey = column[String]("JOB_DETRITUS_KEY")
    def jobDetritusValue = column[String]("JOB_DETRITUS_VALUE")
    def resultMetaInfoId = column[Int]("RESULT_METAINFO_ID")

    override def * = (jobDetritusKey, jobDetritusValue, resultMetaInfoId, callCachingJobDetritusId.?) <>
      (CallCachingJobDetritusEntry.tupled, CallCachingJobDetritusEntry.unapply)

    def callCachingJobDetritusUniquenessConstraint =
      index("UK_CALL_CACHING_JOB_DETRITUS", (jobDetritusKey, resultMetaInfoId), unique = true)

    def callCachingResultMetaInfo = foreignKey(
      "CCJD_RESULT_METAINFO_ID_FK", resultMetaInfoId, callCachingResultMetaInfos)(_.callCachingResultMetaInfoId)
  }

  protected val callCachingJobDetritus = TableQuery[CallCachingJobDetritus]

  val callCachingJobDetritusIdAutoInc =
    callCachingJobDetritus returning callCachingJobDetritus.map(_.callCachingJobDetritusId)

  /**
    * Find all job detrita which match a given RESULT_METAINFO_ID
    */
  val jobDetritusForMetaInfoId = Compiled(
    (resultMetaInfoId: Rep[Int]) => for {
      jobDetritus <- callCachingJobDetritus if jobDetritus.resultMetaInfoId === resultMetaInfoId
    } yield jobDetritus)
}
