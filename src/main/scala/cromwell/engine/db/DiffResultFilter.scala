package cromwell.engine.db

import liquibase.database.Database
import liquibase.diff.{DiffResult, Difference, ObjectDifferences}
import liquibase.structure.DatabaseObject
import liquibase.structure.core.{DataType, Index}

import scala.collection.JavaConverters._

/**
  * Filters liquibase results.
  */
object DiffResultFilter {
  /**
    * A filter for a database object.
    */
  type ObjectFilter = (Database, DatabaseObject) => Boolean

  /**
    * A filter for a database object difference.
    */
  type DiffFilter = (Database, Database, DatabaseObject, Difference) => Boolean

  /**
    * The standard type filters to ignore, including interchangable types.
    */
  val StandardTypeFilters: Seq[DiffFilter] =
    Seq(isTypeSimilar("CLOB", "VARCHAR"), isTypeSimilar("TINYINT", "BOOLEAN"))

  /**
    * Filters diff results using an object filter and changed object filters. Filters that return false are removed.
    *
    * @param diffResult The original diff result.
    * @param changedFilters Filters for changed objects.
    * @param unchangedFilter Filters for missing or unexpected objects.
    * @return The updated diff result.
    */
  def filter(diffResult: DiffResult, changedFilters: Seq[DiffFilter], unchangedFilter: ObjectFilter): DiffResult = {
    val referenceDatabase = diffResult.getReferenceSnapshot.getDatabase
    val comparisonDatabase = diffResult.getComparisonSnapshot.getDatabase

    val missingObjects = diffResult.getMissingObjects.asScala
    val unexpectedObjects = diffResult.getUnexpectedObjects.asScala
    val changedObjects = diffResult.getChangedObjects.asScala

    val newDiffResult = new DiffResult(
      diffResult.getReferenceSnapshot,
      diffResult.getComparisonSnapshot,
      diffResult.getCompareControl)

    missingObjects.filterNot(unchangedFilter(referenceDatabase, _)).foreach(newDiffResult.addMissingObject)
    unexpectedObjects.filterNot(unchangedFilter(comparisonDatabase, _)).foreach(newDiffResult.addUnexpectedObject)

    val filteredChangedObjects =
      changedObjects.filterNot(isSameObject(referenceDatabase, comparisonDatabase, changedFilters))
    for ((obj, difference) <- filteredChangedObjects) {
      newDiffResult.addChangedObject(obj, difference)
    }

    newDiffResult
  }

  /**
    * Given a set of differences as determined by liquibase, runs a set of filters across the set, removing unwanted
    * differences. Once all unwanted differences are removed, if the remaining set is empty, returns true.
    *
    * NOTE: Each filter should return true if the difference should be removed from the set.
    *
    * @param referenceDatabase The reference database.
    * @param comparisonDatabase The comparison database.
    * @param filters Filters to run on each diff.
    * @param objectAndDiff A tuple of the object and the set of differences.
    * @return True if the object is actully the same.
    */
  def isSameObject(referenceDatabase: Database, comparisonDatabase: Database, filters: Seq[DiffFilter])
                  (objectAndDiff: (DatabaseObject, ObjectDifferences)): Boolean = {
    val (obj, objectDifferences) = objectAndDiff
    val differences = objectDifferences.getDifferences.asScala
    val filtered = filters.foldLeft(differences)((diffs, diffFilter) =>
      diffs.filterNot(diffFilter(referenceDatabase, comparisonDatabase, obj, _)))
    filtered.isEmpty
  }

  /**
    * Checks if the difference is due to the schemas being created with different but still similar types.
    *
    * Returns true if the type names are both found in the list, or if the types are the same name.
    *
    * @param similarTypes Types that are similar.
    * @param referenceDatabase The reference database.
    * @param comparisonDatabase The comparison database.
    * @param databaseObject The database object.
    * @param difference The difference reported.
    * @return True if the object is actually similar based on type.
    */
  def isTypeSimilar(similarTypes: String*)
                   (referenceDatabase: Database, comparisonDatabase: Database,
                    databaseObject: DatabaseObject, difference: Difference): Boolean = {
    val compared = difference.getComparedValue
    val referenced = difference.getReferenceValue
    compared.isInstanceOf[DataType] && referenced.isInstanceOf[DataType] && {
      val comparedType = compared.asInstanceOf[DataType].getTypeName
      val referencedType = referenced.asInstanceOf[DataType].getTypeName
      comparedType == referencedType ||
        (similarTypes.contains(comparedType.toUpperCase) && similarTypes.contains(referencedType.toUpperCase))
    }
  }

  /**
    * Returns true if the object is liquibase database object.
    * @param database The source database.
    * @param databaseObject The database object.
    * @return True if the object is a liquibase database object.
    */
  def isLiquibaseObject(database: Database, databaseObject: DatabaseObject): Boolean = {
    database.isLiquibaseObject(databaseObject)
  }

  /**
    * Adds utility methods to a liquibase diff result.
    *
    * @param diffResult The origin diff result.
    */
  implicit class EnhancedDiffResult(val diffResult: DiffResult) extends AnyVal {
    /**
      * Filters liquibase objects from a diff result.
      *
      * @return The diff result without liquibase objects.
      */
    def filterLiquibaseObjects = filter(diffResult, Seq.empty, isLiquibaseObject)

    /**
      * Filters changed objects. Filters that return false are removed.
      *
      * @param filters The filters to apply.
      * @return The updated diff result.
      */
    def filterChangedObjects(filters: Seq[DiffFilter]) = filter(diffResult, filters, (_, _) => false)
  }
}
