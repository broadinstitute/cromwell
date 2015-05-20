package cromwell.server

// Note that as per the language specification, this is instiated lazily and only used when necessary (i.e. server mode)
object CromwellServer extends WorkflowManagerSystem{
  // FIXME: Spray stuff would go here

  println("Cromwell server started")
}

