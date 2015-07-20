package cromwell.engine

// FIXME: Why did I have to move this out of package object?
// FIXME: Was getting error "constructor definition not allowed here"
// FIXME: I think there are interlocking imports on that package object or something like that
//object ExecutionStatus extends Enumeration {
//  type ExecutionStatus = Value
//  val NotStarted, Starting, Running, Failed, Done = Value
//}