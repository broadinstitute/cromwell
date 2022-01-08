package cromwell

object CromwellApp extends App {

  sealed trait Command
  case object Run extends Command
  case object Server extends Command
  case object Submit extends Command
  
  def buildParser(): scopt.OptionParser[CommandLineArguments] = new CommandLineParser()

  def runCromwell(args: CommandLineArguments): Unit = {
    args.command match {
      case Some(Run) => CromwellEntryPoint.runSingle(args)
      case Some(Server) => CromwellEntryPoint.runServer()
      case Some(Submit) => CromwellEntryPoint.submitToServer(args)
      case None => showUsageAndExitWithError()
    }
  }

  val parser = buildParser()

  val parsedArgs = parser.parse(args, CommandLineArguments())
  parsedArgs match {
    case Some(pa) => runCromwell(pa)
    case None => showUsageAndExitWithError()
  }

  private def showUsageAndExitWithError(): Unit = {
    System.err.println(parser.usage)
    System.exit(1)
  }
}
