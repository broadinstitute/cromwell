package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingResultMetaInfoEntry
import slick.profile.RelationalProfile.ColumnOption.Default

import scalaz._

trait CallCachingResultMetaInfoComponent {

  this: DriverComponent with CallCachingHashComponent =>

  import driver.api._

  class CallCachingResultMetaInfoEntries(tag: Tag) extends Table[CallCachingResultMetaInfoEntry](tag, "CALL_CACHING_RESULT_METAINFO") {
    def callCachingResultMetaInfoId = column[Int]("CALL_CACHING_RESULT_METAINFO_ID", O.PrimaryKey, O.AutoInc)

    def workflowUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def callFqn = column[String]("CALL_FQN")

    def returnCode = column[Option[Int]]("RETURN_CODE")

    def scatterIndex = column[Int]("JOB_SCATTER_INDEX")

    def allowResultReuse = column[Boolean]("ALLOW_RESULT_REUSE", Default(true))

    override def * = (workflowUuid, callFqn, scatterIndex, returnCode, allowResultReuse, callCachingResultMetaInfoId.?) <>
      (CallCachingResultMetaInfoEntry.tupled, CallCachingResultMetaInfoEntry.unapply)

    def ccmUniquenessConstraint = index("UK_CALL_CACHING_RESULT_METAINFO", (workflowUuid, callFqn, scatterIndex), unique = true)
  }

  protected val callCachingResultMetaInfos = TableQuery[CallCachingResultMetaInfoEntries]

  val callCachingResultMetaInfoIdsAutoInc =
    callCachingResultMetaInfos returning callCachingResultMetaInfos.map(_.callCachingResultMetaInfoId)

  /**
    * Useful for finding the call caching result meta info for a given ID
    */
  val metaInfoById = Compiled(
    (metaInfoId: Rep[Int]) => for {
      metaInfo <- callCachingResultMetaInfos
      if metaInfo.callCachingResultMetaInfoId === metaInfoId
    } yield metaInfo)

  /**
    * Returns true if there exists a row in callCachingHashes that matches the parameters.
    *
    * @param callCachingResultMetaInfoId The foreign key.
    * @param hashKeyHashValue            The hash key and hash value as a tuple.
    * @return True or false if the row exists.
    */
  private def existsCallCachingResultMetaInfoIdHashKeyHashValue(callCachingResultMetaInfoId: Rep[Int])
                                                               (hashKeyHashValue: (String, String)): Rep[Boolean] = {
    callCachingHashes.filter(callCachingHashEntry =>
      callCachingHashEntry.resultMetaInfoId === callCachingResultMetaInfoId &&
        callCachingHashEntry.hashKey === hashKeyHashValue._1 &&
        callCachingHashEntry.hashValue === hashKeyHashValue._2
    ).exists
  }

  /**
    * Returns true if there exists rows for all the hashKeyHashValues.
    *
    * @param callCachingResultMetaInfoId The foreign key.
    * @param hashKeyHashValues           The hash keys and hash values as tuples.
    * @return True or false if all the rows exist.
    */
  private def existsAllCallCachingResultMetaInfoIdHashKeyHashValues(callCachingResultMetaInfoId: Rep[Int],
                                                                    hashKeyHashValues: NonEmptyList[(String, String)]):
  Rep[Boolean] = {
    hashKeyHashValues.
      map(existsCallCachingResultMetaInfoIdHashKeyHashValue(callCachingResultMetaInfoId)).
      list.toList.reduce(_ && _)
  }

  /**
    * Returns the callCachingResultMetaInfoId where all the hash keys and hash values exist.
    *
    * @param hashKeyHashValues The hash keys and hash values as tuples.
    * @return The callCachingResultMetaInfoIds with all of the hash keys and hash values.
    */
  def callCachingResultMetaInfoIdByHashKeyHashValues(hashKeyHashValues: NonEmptyList[(String, String)]) = {
    for {
      callCachingResultMetaInfo <- callCachingResultMetaInfos
      if existsAllCallCachingResultMetaInfoIdHashKeyHashValues(
        callCachingResultMetaInfo.callCachingResultMetaInfoId, hashKeyHashValues)
    } yield callCachingResultMetaInfo.callCachingResultMetaInfoId
  }
}
