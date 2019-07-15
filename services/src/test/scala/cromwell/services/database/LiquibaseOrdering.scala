package cromwell.services.database

import liquibase.structure.core._

object LiquibaseOrdering {
  implicit val liquibaseOrderingRelation: Ordering[Relation] = Ordering.by[Relation, String](_.getName)

  implicit val liquibaseOrderingTable: Ordering[Table] = Ordering.by[Table, Relation](identity)

  implicit val liquibaseOrderingColumn: Ordering[Column] = Ordering.by[Column, String] { column =>
    s"${column.getRelation.getName}.${column.getName}"
  }

  implicit val liquibaseOrderingPrimaryKey: Ordering[PrimaryKey] = Ordering.by[PrimaryKey, String](_.getName)

  implicit val liquibaseOrderingForeignKey: Ordering[ForeignKey] = Ordering.by[ForeignKey, String](_.getName)

  implicit val liquibaseOrderingUniqueConstraint: Ordering[UniqueConstraint] = {
    Ordering.by[UniqueConstraint, String](_.getName)
  }

  implicit val liquibaseOrderingIndex: Ordering[Index] = Ordering.by[Index, String](_.getName)
}
