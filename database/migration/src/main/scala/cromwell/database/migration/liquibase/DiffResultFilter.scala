package cromwell.database.migration.liquibase

import liquibase.database.Database
import liquibase.diff.{DiffResult, Difference, ObjectDifferences}
import liquibase.structure.DatabaseObject
import liquibase.structure.core._

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
    * The standard type filters to ignore, including interchangeable types.
    */
  val StandardTypeFilters: Seq[DiffFilter] = Seq(isReordered, isVarchar255, isTypeSimilar("TINYINT", "BOOLEAN"))

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
    * @return True if the object is actually the same.
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
    * Checks if the difference is due to the schemas being created with our default varchar(255).
    *
    * @param referenceDatabase  The reference database.
    * @param comparisonDatabase The comparison database.
    * @param databaseObject     The database object.
    * @param difference         The difference reported.
    * @return True if the object is actually the same with slightly different column widths.
    */
  def isVarchar255(referenceDatabase: Database, comparisonDatabase: Database,
                   databaseObject: DatabaseObject, difference: Difference): Boolean = {
    val compared = difference.getComparedValue
    val referenced = difference.getReferenceValue
    compared.isInstanceOf[DataType] && referenced.isInstanceOf[DataType] && {
      val comparedDataType = compared.asInstanceOf[DataType]
      val referencedDataType = referenced.asInstanceOf[DataType]
      comparedDataType.getTypeName == "VARCHAR" && referencedDataType.getTypeName == "VARCHAR" &&
        // Our liquibase copypasta defaults VARCHAR to 255. Slick without a value defaults to 254
        (comparedDataType.getColumnSize + referencedDataType.getColumnSize == 255 + 254)
    }
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
      val comparedType = compared.asInstanceOf[DataType].getTypeName.toUpperCase
      val referencedType = referenced.asInstanceOf[DataType].getTypeName.toUpperCase
      similarTypes.contains(comparedType) && similarTypes.contains(referencedType.toUpperCase)
    }
  }

  /**
    * Checks if the difference is due to the schemas being reordered.
    *
    * Returns true if the type is reordered.
    *
    * @param referenceDatabase The reference database.
    * @param comparisonDatabase The comparison database.
    * @param databaseObject The database object.
    * @param difference The difference reported.
    * @return True if the object is actually similar based on type.
    */
  def isReordered(referenceDatabase: Database, comparisonDatabase: Database,
                  databaseObject: DatabaseObject, difference: Difference): Boolean = {
    difference.getField == "order"
  }

  /**
    * Returns true if the object is a change log object.
    *
    * @param database The source database.
    * @param databaseObject The database object.
    * @return True if the object is a change log object.
    */
  def isChangeLog(database: Database, databaseObject: DatabaseObject): Boolean = {
    databaseObject match {
      case table: Table => table.getName.contains("DATABASECHANGELOG")
      case column: Column => isChangeLog(database, column.getRelation)
      case index: Index => isChangeLog(database, index.getRelation)
      case key: PrimaryKey => isChangeLog(database, key.getTable)
      case _ => false
    }
  }

  /**
    * Returns true if the object is liquibase database object.
    *
    * @param database The source database.
    * @param databaseObject The database object.
    * @return True if the object is a liquibase database object.
    */
  def isLiquibaseObject(database: Database, databaseObject: DatabaseObject): Boolean = {
    database.isLiquibaseObject(databaseObject)
  }

  /**
    * Returns true if the object is a member of the excluded table.
    *
    * @param tables         Tables to check.
    * @param database       The source database.
    * @param databaseObject The database object.
    * @return True if the object is a member of the tables.
    */
  def isTableObject(tables: Seq[String])
                   (database: Database, databaseObject: DatabaseObject): Boolean = {
    isTableObject(tables, databaseObject)
  }

  private def isTableObject(tables: Seq[String], databaseObject: DatabaseObject): Boolean = {
    tables.exists(table =>
      databaseObject.getName.equalsIgnoreCase(table) ||
        getContainingObjects(databaseObject).exists(isTableObject(tables, _))
    )
  }

  // getContainingObjects is ill-mannered and returns null when really it ought to return an empty array, so wrap
  // in an `Option` and `getOrElse`.
  private def getContainingObjects(databaseObject: DatabaseObject): Array[DatabaseObject] = {
    Option(databaseObject.getContainingObjects).getOrElse(Array.empty)
  }

  /**
    * Adds utility methods to a liquibase diff result.
    *
    * @param diffResult The origin diff result.
    */
  implicit class EnhancedDiffResult(val diffResult: DiffResult) extends AnyVal {
    /**
      * Filters changelogs.
      *
      * @return The diff result without changelogs.
      */
    def filterChangeLogs = filter(diffResult, Seq.empty, isChangeLog)

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

    def filterTableObjects(tables: Seq[String]) = filter(diffResult, Seq.empty, isTableObject(tables))
  }
}
