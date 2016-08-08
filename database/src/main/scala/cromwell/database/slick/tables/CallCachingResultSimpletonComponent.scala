package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingResultSimpletonEntry

trait CallCachingResultSimpletonComponent {

  this: DriverComponent with CallCachingResultMetaInfoComponent =>

  import driver.api._

  class CallCachingResultSimpletonEntries(tag: Tag) extends Table[CallCachingResultSimpletonEntry](tag, "CALL_CACHING_RESULT_SIMPLETON") {
    def callCachingResultSimpletonId = column[Int]("CALL_CACHING_RESULT_SIMPLETON_ID", O.PrimaryKey, O.AutoInc)
    def simpletonKey = column[String]("SIMPLETON_KEY")
    def simpletonValue = column[String]("SIMPLETON_VALUE")
    def wdlType = column[String]("WDL_TYPE")
    def resultMetaInfoId = column[Int]("RESULT_METAINFO_ID")

    override def * = (simpletonKey, simpletonValue, wdlType, resultMetaInfoId, callCachingResultSimpletonId.?) <>
      (CallCachingResultSimpletonEntry.tupled, CallCachingResultSimpletonEntry.unapply)

    def ccrsUniquenessConstraint = index("UK_CALL_CACHING_RESULT_SIMPLETON", (simpletonKey, resultMetaInfoId), unique = true)
    def resultMetaInfoForeignKey = foreignKey("CCRS_RESULT_METAINFO_ID_FK", resultMetaInfoId, callCachingResultMetaInfo)(_.callCachingResultMetaInfoId)
  }

  protected val callCachingResultSimpletons = TableQuery[CallCachingResultSimpletonEntries]

  val callCachingResultSimpletonAutoInc = callCachingResultSimpletons returning callCachingResultSimpletons.map(_.callCachingResultSimpletonId)

  /**
    * Find all result simpletons which match a given RESULT_METAINFO_ID
    */
  val resultSimpletonsForMetaInfoId = Compiled(
    (resultMetaInfoId: Rep[Int]) => for {
      simpleton <- callCachingResultSimpletons if simpleton.resultMetaInfoId === resultMetaInfoId
    } yield simpleton)
}
