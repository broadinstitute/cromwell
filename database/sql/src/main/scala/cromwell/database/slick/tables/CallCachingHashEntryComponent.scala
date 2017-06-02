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

    def ucCallCachingHashEntryCceiHk =
      index("UC_CALL_CACHING_HASH_ENTRY_CCEI_HK", (callCachingEntryId, hashKey), unique = true)
  }

  val callCachingHashEntries = TableQuery[CallCachingHashEntries]

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
    * Returns whether or not there exists at least one callCachingEntryId for which all the hash keys and hash values match, and are allowed to be reused.
    */
  def existsMatchingCachingEntryIdsForHashKeyHashValues(hashKeyHashValues: NonEmptyList[(String, String)]) = {
    (for {
      callCachingEntry <- callCachingEntries
      if existsAllCallCachingEntryIdHashKeyHashValues(
        callCachingEntry.callCachingEntryId, hashKeyHashValues) && callCachingEntry.allowResultReuse
    } yield callCachingEntry.callCachingEntryId).exists
  }

  val callCachingEntriesForWorkflowFqnIndex = Compiled(
    (
      callA: (Rep[String], Rep[String],  Rep[Int]),
      callB: (Rep[String], Rep[String], Rep[Int])
    ) => {

      // Query that returns the hashes associated with a workflowExecutionUuid / callFqn / jobIndex triplet
      // The attempt is not part of the filter, which means if there are several attempts stored, all hashes will be returned.
      // At the time where this is written, hashes for failed (but retried) attempts are not stored. If that behavior were
      // to change this logic would probably need to be updated.
      def makeHashQuery(call: (Rep[String], Rep[String],  Rep[Int])) = {
        val (workflowExecutionUuid, callFqn, jobIndex) = call
        for {
          callCachingEntry <- callCachingEntries
          if callCachingEntry.workflowExecutionUuid === workflowExecutionUuid
          if callCachingEntry.callFullyQualifiedName === callFqn
          if callCachingEntry.jobIndex === jobIndex
          callCacheHashes <- callCachingHashEntries
          if callCacheHashes.callCachingEntryId === callCachingEntry.callCachingEntryId
        } yield callCacheHashes
      }

      // Hashes for call A
      val hashEntriesForA = makeHashQuery(callA)
      // Hashes for call B
      val hashEntriesForB = makeHashQuery(callB)

      for {
        hashes <- hashEntriesForA
          // Join both hash sets. Full join so we get everything, including hashes in A but not in B and vice versa
          .joinFull(hashEntriesForB)
          // Join on hashKeys
          .on(_.hashKey === _.hashKey)
          // Only keep hashes for which the values are not the same
          // Note that because we're dealing with Rep[Option[...]] and not Option[...] 
          // we can't pattern match the maybeHashes to Some or None
          .filter({ 
          case (maybeHashA, maybeHashB) =>
            // HashKey exists in B but not in A
            (maybeHashA.isEmpty && maybeHashB.nonEmpty) ||
            // HashKey exists in A but not in B
            (maybeHashA.nonEmpty && maybeHashB.isEmpty) ||
            // HashKey is in both but has different values
            (for {
              hashA <- maybeHashA
              hashB <- maybeHashB
            } yield hashA.hashValue =!= hashB.hashValue)
            // Both A and B are null, which should not be possible...  
            .getOrElse(false)
          })
          // Remove duplicates
          .distinct
          // Project only hashKey -> hashValue pairs
          .map({
            case (maybeHashA, maybeHashB) =>
              maybeHashA.map(s => s.hashKey -> s.hashValue) -> maybeHashB.map(s => s.hashKey -> s.hashValue)
          })
      } yield hashes
    }
  )
}
