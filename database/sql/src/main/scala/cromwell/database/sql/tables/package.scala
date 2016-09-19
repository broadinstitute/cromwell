package cromwell.database.sql

/**
  * The collection of case classes representing database tables.
  *
  * It's possible that we may need to swap out the database layers at some point, because:
  * - slick upgrades from 3.x to 4.x (see previous 2.x to 3.x migration)
  * - we want a new relational mapping layer
  * - we simply want a mock database
  *
  * '''TL;DR Keep the database logic to an absolute minimum!'''
  *
  * Each case class name should match the table name, replacing capital letters in the class name with an underscore,
  * and then converting the entire string to upper case.
  *
  * The columns in each class should be a primitive type matching the database column.
  *
  * All column types should match the type in the database, and be one of:
  * - `Boolean`
  * - `Double`
  * - `Int`
  * - `Long`
  * - `java.sql.Clob`
  * - `java.sql.Timestamp`
  *
  * Nullable columns should be wrapped in an `Option`.
  *
  * Primary and foreign key columns are the only columns that should be defaulted, as they are to be filled in by the
  * database, and cannot and should not be set within the business logic. On the other hand, columns to be filled in by
  * the business logic should __not__ have a default in the database package, even if they are nullable.
  *
  * Example:
  *
  * {{{
  * case class Car
  * (
  *   make: String,              // Generic make as a String. No enums, traits, objects, or other business logic.
  *   model: String,             // Same for model. Any tracking of "make must belong to model" isn't done in ORM/FRM.
  *   year: Int,                 // Generic year here. Any business logic handled elsewhere.
  *   lastOwner: Option[String]  // No defaults! Please let the business layer populate this value.
  *   carId: Option[Int] = None  // PK, will be automatically populated by the database. Last because it's defaulted.
  * )
  * }}}
  *
  * The database(s) will store whatever model passed in. Place the common car types in core, not in the database layer.
  *
  * {{{
  * enum CarModels { Ford, Toyota, Volkswagen, ... } <<-- NOT IN THE DATABASE!
  * }}}
  *
  * Then, to generate a car in the business layer, all values to be passed into the database must be specified,
  * including an explicit statement about the last owner:
  *
  * {{{
  * val newCar = Car(Ford.toString, Ford.Fusion.toString, 2010, None)
  * }}}
  */
package object tables
