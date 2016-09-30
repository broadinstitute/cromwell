package cromwell.database.slick.tables

import cats.data.NonEmptyList
import cromwell.database.sql.tables.CallCachingHashEntry

trait CallCachingHashEntryComponent {

  this: DriverComponent with CallCachingEntryComponent =>

  import driver.api._

  class CallCachingHashEntries(tag: Tag) extends Table[CallCachingHashEntry](tag, "CALL_CACHING_HASH_ENTRY") {
    def callCachingHashEntryId = column[Int]("CALL_CACHING_HASH_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def hashKey = column[String]("HASH_KEY")

    def hashValue = column[String]("HASH_VALUE")

    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID")

    override def * = (hashKey, hashValue, callCachingEntryId.?, callCachingHashEntryId.?) <>
      (CallCachingHashEntry.tupled, CallCachingHashEntry.unapply)

    def fkCallCachingHashEntryCallCachingEntryId = foreignKey("FK_CALL_CACHING_HASH_ENTRY_CALL_CACHING_ENTRY_ID",
      callCachingEntryId, callCachingEntries)(_.callCachingEntryId)

    def ucCallCachingHashEntryCcei =
      index("UC_CALL_CACHING_HASH_ENTRY_CCEI", (callCachingEntryId, hashKey), unique = true)
  }

  protected val callCachingHashEntries = TableQuery[CallCachingHashEntries]

  val callCachingHashEntryIdsAutoInc = callCachingHashEntries returning
    callCachingHashEntries.map(_.callCachingHashEntryId)

  /**
    * Returns true if there exists a row in callCachingHashEntries that matches the parameters.
    *
    * @param callCachingEntryId The foreign key.
    * @param hashKeyHashValue   The hash key and hash value as a tuple.
    * @return True or false if the row exists.
    */
  private def existsCallCachingEntryIdHashKeyHashValue(callCachingEntryId: Rep[Int])
                                                      (hashKeyHashValue: (String, String)): Rep[Boolean] = {
    callCachingHashEntries.filter(callCachingHashEntry =>
      callCachingHashEntry.callCachingEntryId === callCachingEntryId &&
        callCachingHashEntry.hashKey === hashKeyHashValue._1 &&
        callCachingHashEntry.hashValue === hashKeyHashValue._2
    ).exists
  }

  /**
    * Returns true if there exists rows for all the hashKeyHashValues.
    *
    * @param callCachingEntryId The foreign key.
    * @param hashKeyHashValues  The hash keys and hash values as tuples.
    * @return True or false if all the rows exist.
    */
  private def existsAllCallCachingEntryIdHashKeyHashValues(callCachingEntryId: Rep[Int],
                                                           hashKeyHashValues: NonEmptyList[(String, String)]):
  Rep[Boolean] = {
    hashKeyHashValues.
      map(existsCallCachingEntryIdHashKeyHashValue(callCachingEntryId)).
      toList.reduce(_ && _)
  }

  /**
    * Returns the callCachingEntryIds where all the hash keys and hash values exist.
    *
    * @param hashKeyHashValues The hash keys and hash values as tuples.
    * @return The callCachingEntryIds with all of the hash keys and hash values.
    */
  def callCachingEntryIdsForHashKeyHashValues(hashKeyHashValues: NonEmptyList[(String, String)]) = {
    for {
      callCachingEntry <- callCachingEntries
      if existsAllCallCachingEntryIdHashKeyHashValues(
        callCachingEntry.callCachingEntryId, hashKeyHashValues)
    } yield callCachingEntry.callCachingEntryId
  }
}
