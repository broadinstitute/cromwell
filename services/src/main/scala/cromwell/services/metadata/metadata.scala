package cromwell.services.metadata

case class QueryParameter(key: String, value: String)

object Patterns {
  val WorkflowName = """
        (?x)                             # Turn on comments and whitespace insensitivity.

        (                                # Begin capture.

          [a-zA-Z][a-zA-Z0-9_]*          # WDL identifier naming pattern of an initial alpha character followed by zero
                                         # or more alphanumeric or underscore characters.

        )                                # End capture.
                     """.trim.r

  val CallFullyQualifiedName = """
      (?x)                               # Turn on comments and whitespace insensitivity.

      (                                  # Begin outer capturing group for FQN.

        (?:[a-zA-Z][a-zA-Z0-9_]*)        #   Inner noncapturing group for top-level workflow name. This is the WDL
                                         #   identifier naming pattern of an initial alpha character followed by zero
                                         #   or more alphanumeric or underscore characters.

        (?:\.[a-zA-Z][a-zA-Z0-9_]*){1}   #   Inner noncapturing group for call name, a literal dot followed by a WDL
                                         #   identifier.  Currently this is quantified to {1} since the call name is
                                         #   mandatory and nested workflows are not supported.  This could be changed
                                         #   to + or a different quantifier if these assumptions change.

      )                                  # End outer capturing group for FQN.


      (?:                                # Begin outer noncapturing group for shard.
        \.                               #   Literal dot.
        (\d+)                            #   Captured shard digits.
      )?                                 # End outer optional noncapturing group for shard.
                               """.trim.r                         // The trim is necessary as (?x) must be at the beginning of the regex.
}
