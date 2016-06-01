package cromwell

import java.nio.file.Paths

import akka.actor.{Actor, ActorRef, ActorSystem, Props}
import akka.pattern.ask
import akka.testkit._
import akka.util.Timeout
import better.files.File
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestkitSpec._
import cromwell.backend._
import cromwell.core.WorkflowId
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine._
import cromwell.engine.backend.{BackendConfigurationEntry, CallLogs}
import cromwell.engine.workflow.{WorkflowManagerActor, WorkflowMetadataKeys}
import cromwell.server.WorkflowManagerSystem
import cromwell.services.MetadataServiceActor._
import cromwell.services.{MetadataQuery, ServiceRegistryClient}
import cromwell.util.SampleWdl
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.MetadataBuilderActor
import cromwell.webservice.PerRequest.RequestComplete
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, Matchers, OneInstancePerTest, WordSpecLike}
import spray.json._
import wdl4s.Call
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}
import wdl4s.types.{WdlArrayType, WdlMapType, WdlStringType}
import wdl4s.values._

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.language.postfixOps
import scala.reflect.ClassTag
import scala.util.matching.Regex

case class TestBackendLifecycleActorFactory(config: Config) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                configurationDescriptor: BackendConfigurationDescriptor): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      configurationDescriptor: BackendConfigurationDescriptor): Props = {
    throw new NotImplementedError("this is not implemented")
  }

  override def workflowFinalizationActorProps(): Option[Props] = None

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           configurationDescriptor: BackendConfigurationDescriptor): WdlStandardLibraryFunctions = {
    NoFunctions
  }
}

case class OutputNotFoundException(outputFqn: String, actualOutputs: String) extends RuntimeException(s"Expected output $outputFqn was not found in: '$actualOutputs'")

object CromwellTestkitSpec {
  val ConfigText =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "INFO"
      |  actor {
      |    debug {
      |       receive = on
      |    }
      |  }
      |  dispatchers {
      |    slow-actor-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |  }
      |  test {
      |    # Some of our tests fire off a message, then expect a particular event message within 3s (the default).
      |    # Especially on CI, the metadata test does not seem to be returning in time. So, overriding the timeouts
      |    # with slightly higher values. Alternatively, could also adjust the akka.test.timefactor only in CI.
      |    filter-leeway = 5s
      |    single-expect-default = 5s
      |    default-timeout = 10s
      |  }
      |}
    """.stripMargin

  val timeoutDuration = 30 seconds

  class TestWorkflowManagerSystem extends WorkflowManagerSystem {
    override protected def systemName: String = "test-system"
    override protected def newActorSystem() = ActorSystem(systemName, ConfigFactory.parseString(CromwellTestkitSpec.ConfigText))
    /**
      * Do NOT shut down the test actor system inside the normal flow.
      * The actor system will be externally shutdown outside the block.
      */
    override def shutdownActorSystem() = {}

    def shutdownTestActorSystem() = super.shutdownActorSystem()
  }

  /**
    * Loans a test actor system. NOTE: This should be run OUTSIDE of a wait block, never within one.
    */
  def withTestWorkflowManagerSystem[T](block: WorkflowManagerSystem => T): T = {
    val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem
    try {
      block(testWorkflowManagerSystem)
    } finally {
      TestKit.shutdownActorSystem(testWorkflowManagerSystem.actorSystem, timeoutDuration)
    }
  }

  /**
   * Wait for exactly one occurrence of the specified info pattern in the specified block.  The block is in its own
   * parameter list for usage syntax reasons.
   */
  def waitForInfo[T](pattern: String, occurrences: Int = 1)(block: => T)(implicit system: ActorSystem): T = {
    EventFilter.info(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  /**
   * Wait for occurrence(s) of the specified warning pattern in the specified block.  The block is in its own parameter
   * list for usage syntax reasons.
   */
  def waitForWarning[T](pattern: String, occurrences: Int = 1)(block: => T)(implicit system: ActorSystem): T = {
    EventFilter.warning(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  /**
   * Wait for occurrence(s) of the specified error pattern in the specified block.  The block is in its own parameter
   * list for usage syntax reasons.
   */
  def waitForError[T](pattern: String, occurrences: Int = 1)(block: => T)(implicit system: ActorSystem): T = {
    EventFilter.error(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  /**
   * Akka TestKit appears to be unable to match errors generated by `log.error(Throwable, String)` with the normal
   * `EventFilter.error(...).intercept {...}` mechanism since `EventFilter.error` forces the use of a dummy exception
   * that never matches a real exception.  This method works around that problem by building an `ErrorFilter` more
   * explicitly to allow the caller to specify a `Throwable` class.
   */
  def waitForErrorWithException[T](pattern: String,
                                   throwableClass: Class[_ <: Throwable] = classOf[Throwable],
                                   occurrences: Int = 1)
                                  (block: => T)
                                  (implicit system: ActorSystem): T = {
    val regex = Right[String, Regex](pattern.r)
    ErrorFilter(throwableClass, source = None, message = regex, complete = false)(occurrences = occurrences).intercept {
      block
    }
  }

  /**
    * Special case for validating outputs. Used when the test wants to check that an output exists, but doesn't care what
    * the actual value was.
    */
  lazy val AnyValueIsFine: WdlValue = WdlString("Today you are you! That is truer than true! There is no one alive who is you-er than you!")

  implicit class EnhancedWorkflowManagerActor(val manager: TestActorRef[WorkflowManagerActor]) extends AnyVal {
    // `implicit` for the asks below.
    private implicit def timeout: Timeout = 30 seconds

    def submit(sources: WorkflowSourceFiles): WorkflowId = {
      val submitMessage = WorkflowManagerActor.SubmitWorkflowCommand(sources)
      Await.result(manager.ask(submitMessage), Duration.Inf).asInstanceOf[WorkflowManagerSubmitSuccess].id
    }

    @deprecated("callers should use getWorkflowOutputs instead, but the return type has changed")
    def workflowOutputs(id: WorkflowId): engine.WorkflowOutputs = throw new UnsupportedOperationException()

    @deprecated("should no longer be on WorkflowManagerActor in a PBE world")
    def workflowStdoutStderr(id: WorkflowId): Map[FullyQualifiedName, Seq[CallLogs]] = throw new UnsupportedOperationException()

    @deprecated("should no longer be on WorkflowManagerActor in a PBE world")
    def callStdoutStderr(id: WorkflowId, fqn: FullyQualifiedName): Seq[CallLogs] = throw new UnsupportedOperationException()
  }

  lazy val DefaultConfig = ConfigFactory.load
  lazy val DefaultLocalBackendConfig = ConfigFactory.parseString(
    """
      |  {
      |    root: "cromwell-executions"
      |
      |    filesystems {
      |      local {
      |        localization: [
      |          "hard-link", "soft-link", "copy"
      |        ]
      |      }
      |    }
      |  }
    """.stripMargin)

  lazy val DefaultLocalBackendConfigEntry = BackendConfigurationEntry(
    name = "local",
    lifecycleActorFactoryClass = "cromwell.TestBackendLifecycleActorFactory",
    DefaultLocalBackendConfig
  )

  lazy val JesBackendConfig = ConfigFactory.parseString(
    """
      |{
      |  // Google project
      |  project = "my-cromwell-workflows"
      |
      |  // Base bucket for workflow executions
      |  root = "gs://my-cromwell-workflows-bucket"
      |
      |  // Polling for completion backs-off gradually for slower-running jobs.
      |  // This is the maximum polling interval (in seconds):
      |  maximum-polling-interval = 600
      |
      |  // Optional Dockerhub Credentials. Can be used to access private docker images.
      |  dockerhub {
      |    // account = ""
      |    // token = ""
      |  }
      |
      |  genomics {
      |    // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
      |    // Pipelines and manipulate auth JSONs.
      |    auth = "service-account"
      |
      |    // Endpoint for APIs, no reason to change this unless directed by Google.
      |    endpoint-url = "https://genomics.googleapis.com/"
      |  }
      |
      |  filesystems {
      |    gcs {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "user-via-refresh"
      |    }
      |  }
      |}
    """.stripMargin)

  lazy val JesBackendConfigEntry = BackendConfigurationEntry(
    name = "JES",
    lifecycleActorFactoryClass = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory",
    JesBackendConfig
  )
}

abstract class CromwellTestkitSpec extends TestKit(new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)
  with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll with ScalaFutures with OneInstancePerTest with ServiceRegistryClient {

  val name = this.getClass.getSimpleName
  implicit val defaultPatience = PatienceConfig(timeout = Span(30, Seconds), interval = Span(100, Millis))
  implicit val ec = system.dispatcher

  def startingCallsFilter[T](callNames: String*)(block: => T): T =
    waitForPattern(s"starting calls: ${callNames.mkString(", ")}$$") {
      block
    }

  def waitForHandledMessage[T](named: String)(block: => T): T = {
    waitForHandledMessagePattern(s"^received handled message $named") {
      block
    }
  }

  def waitForHandledMessagePattern[T](pattern: String)(block: => T): T = {
    waitForInfo(pattern = pattern, occurrences = 1) {
      block
    }
  }

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Downcasts the `Future[Any]` to a `Future[M]`.
   * <li> Issues a blocking `Await.result` on the `Future`, yielding an `M`.
   * </ol>
   *
   */
  def messageAndWait[M: ClassTag](message: AnyRef)(implicit actorRef: ActorRef): M = {
    val futureAny = actorRef ? message
    Await.result(futureAny.mapTo[M], timeoutDuration)
  }

  /**
   * Wait for exactly one occurrence of the specified info pattern in the specified block within a limited amount of
   * time. The block is in its own parameter list for usage syntax reasons.
   */
  def waitForPattern[T](pattern: String, occurrences: Int = 1)(block: => T): T = {
    within(timeoutDuration) {
      waitForInfo(pattern, occurrences) {
        block
      }
    }
  }

  private def buildWorkflowManagerActor(config: Config) = {
    TestActorRef(new WorkflowManagerActor(config))
  }

  // output received from the metadata service is untyped, so WdlValues are converted to strings
  private def validateOutput(output: WdlValue, expectedOutput: WdlValue): Unit = expectedOutput match {
    case expectedFile: WdlFile => Paths.get(output.valueString).getFileName shouldEqual Paths.get(expectedOutput.valueString).getFileName
    case expectedArray: WdlArray if output.isInstanceOf[WdlArray] =>
      val actualArray = output.asInstanceOf[WdlArray]
      actualArray.value.size should be(expectedArray.value.size)
      (actualArray.value zip expectedArray.value) foreach {
        case (actual, expected) => validateOutput(actual, expected)
      }
    case _ =>
      output.wdlType shouldEqual expectedOutput.wdlType
      output.valueString shouldEqual expectedOutput.valueString
  }

  def runWdl(sampleWdl: SampleWdl,
             eventFilter: EventFilter,
             runtime: String = "",
             workflowOptions: String = "{}",
             terminalState: WorkflowState = WorkflowSucceeded,
             config: Config = DefaultConfig)(implicit ec: ExecutionContext): Map[FullyQualifiedName, WdlValue] = {
    val wma = buildWorkflowManagerActor(config)
    val sources = WorkflowSourceFiles(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, workflowOptions)
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = wma.submit(sources)
        verifyWorkflowState(wma, workflowId, terminalState)
        getWorkflowOutputsFromMetadata(workflowId)
      }
    }
  }

  def runWdlAndAssertOutputs(sampleWdl: SampleWdl,
                             eventFilter: EventFilter,
                             expectedOutputs: Map[FullyQualifiedName, WdlValue],
                             runtime: String = "",
                             workflowOptions: String = "{}",
                             allowOtherOutputs: Boolean = true,
                             terminalState: WorkflowState = WorkflowSucceeded,
                             config: Config = DefaultConfig,
                             workflowManagerActor: Option[TestActorRef[WorkflowManagerActor]] = None)
                            (implicit ec: ExecutionContext): WorkflowId = {
    val workflowManager = workflowManagerActor.getOrElse(buildWorkflowManagerActor(config))
    val sources = sampleWdl.asWorkflowSources(runtime, workflowOptions)
    val maxRetries = 3
    def isFatal(e: Throwable) = e match {
      case _: OutputNotFoundException => false
      case _ => true
    }
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = workflowManager.submit(sources)
        verifyWorkflowState(workflowManager, workflowId, terminalState)

        def checkOutputs() = Future {
          val outputs = getWorkflowOutputsFromMetadata(workflowId)
          val actualOutputNames = outputs.keys mkString ", "
          val expectedOutputNames = expectedOutputs.keys mkString " "

          expectedOutputs foreach { case (outputFqn, expectedValue) =>
            val actualValue = outputs.getOrElse(outputFqn, throw new OutputNotFoundException(outputFqn, actualOutputNames))
            if (expectedValue != AnyValueIsFine) validateOutput(actualValue, expectedValue)
          }
          if (!allowOtherOutputs) {
            outputs foreach { case (actualFqn, actualValue) =>
              val expectedValue = expectedOutputs.getOrElse(actualFqn, throw new RuntimeException(s"Actual output $actualFqn was not wanted in '$expectedOutputNames'"))
              if (expectedValue != AnyValueIsFine) validateOutput(actualValue, expectedValue)
            }
          }
          workflowId
        }

        // Retry because we have no guarantee that all the metadata is there when the workflow finishes...
        Await.result(Retry.withRetry(
          checkOutputs,
          maxRetries = Some(maxRetries),
          isFatal = isFatal,
          backoff = SimpleExponentialBackoff(0.5 seconds, 0.5 second, 1D)), timeoutDuration)
      }
    }
  }

  /*
     FIXME: I renamed this as it appears to be asserting the stdout/stderr of a single call which is kinda weird for
     a full workflow type of thing
  */
  def runSingleCallWdlWithWorkflowManagerActor(wma: TestActorRef[WorkflowManagerActor],
                                               sources: WorkflowSourceFiles,
                                               eventFilter: EventFilter,
                                               fqn: FullyQualifiedName,
                                               index: ExecutionIndex,
                                               stdout: Option[Seq[String]],
                                               stderr: Option[Seq[String]],
                                               expectedOutputs: Map[FullyQualifiedName, WdlValue] = Map.empty)(implicit ec: ExecutionContext) = {
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = wma.submit(sources)
        verifyWorkflowState(wma, workflowId, WorkflowSucceeded)
        val standardStreams = wma.callStdoutStderr(workflowId, fqn)
        stdout foreach { souts =>
          souts shouldEqual (standardStreams map { s => File(s.stdout.value).contentAsString })
        }
        stderr foreach { serrs =>
          serrs shouldEqual (standardStreams map { s => File(s.stderr.value).contentAsString })
        }
      }
    }
  }

  def runWdlWithWorkflowManagerActor(wma: TestActorRef[WorkflowManagerActor],
                                     sources: WorkflowSourceFiles,
                                     eventFilter: EventFilter,
                                     stdout: Map[FullyQualifiedName, Seq[String]],
                                     stderr: Map[FullyQualifiedName, Seq[String]],
                                     expectedOutputs: Map[FullyQualifiedName, WdlValue] = Map.empty,
                                     terminalState: WorkflowState = WorkflowSucceeded,
                                     assertStdoutStderr: Boolean = false)(implicit ec: ExecutionContext) = {
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = wma.submit(sources)
        verifyWorkflowState(wma, workflowId, terminalState)

        if (assertStdoutStderr) {
          val standardStreams = wma.workflowStdoutStderr(workflowId)
          stdout foreach {
            case (fqn, out) if standardStreams.contains(fqn) =>
              out shouldEqual (standardStreams(fqn) map { s => File(s.stdout.value).contentAsString })
          }
          stderr foreach {
            case (fqn, err) if standardStreams.contains(fqn) =>
              err shouldEqual (standardStreams(fqn) map { s => File(s.stderr.value).contentAsString })
          }
        }
      }
    }
  }

  def runWdlAndAssertStdoutStderr(sampleWdl: SampleWdl,
                                  eventFilter: EventFilter,
                                  fqn: FullyQualifiedName,
                                  index: ExecutionIndex,
                                  runtime: String = "",
                                  stdout: Option[Seq[String]] = None,
                                  stderr: Option[Seq[String]] = None,
                                  config: Config = DefaultConfig)
                                 (implicit ec: ExecutionContext) = {
    val actor = buildWorkflowManagerActor(config)
    val sources = WorkflowSourceFiles(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, "{}")
    runSingleCallWdlWithWorkflowManagerActor(actor, sources, eventFilter, fqn, index, stdout, stderr)
  }

  def runWdlAndAssertWorkflowStdoutStderr(sampleWdl: SampleWdl,
                                          eventFilter: EventFilter,
                                          runtime: String = "",
                                          stdout: Map[FullyQualifiedName, Seq[String]] = Map.empty[FullyQualifiedName, Seq[String]],
                                          stderr: Map[FullyQualifiedName, Seq[String]] = Map.empty[FullyQualifiedName, Seq[String]],
                                          terminalState: WorkflowState = WorkflowSucceeded,
                                          config: Config = DefaultConfig)
                                         (implicit ec: ExecutionContext) = {
    val actor = buildWorkflowManagerActor(config)
    val workflowSources = WorkflowSourceFiles(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, "{}")
    runWdlWithWorkflowManagerActor(actor, workflowSources, eventFilter, stdout, stderr, Map.empty, terminalState)
  }

  def getWorkflowMetadata(workflowId: WorkflowId, key: Option[String] = None)(implicit ec: ExecutionContext): JsObject = {
    // MetadataBuilderActor sends it's response to context.parent, so we can't just use an ask to talk to it here
    val supervisor = TestActorRef(new Actor() {
      var originalSender = system.deadLetters

      override def receive: Receive = {
        case m: RequestComplete[JsObject] =>
          originalSender ! m
        case m: GetMetadataQueryAction =>
          originalSender = sender
          val mdba = TestActorRef(MetadataBuilderActor.props(serviceRegistryActor), self, "mdba" )
          mdba ! m
      }
    })

    val message = GetMetadataQueryAction(MetadataQuery(workflowId, None, key))
    Await.result(supervisor.ask(message).mapTo[RequestComplete[(String,JsObject)]], Duration.Inf).response._2
  }

  def verifyWorkflowState(wma: ActorRef, workflowId: WorkflowId, expectedState: WorkflowState)(implicit ec: ExecutionContext): Unit = {
    // Continuously check the state of the workflow until it is in a terminal state
    awaitCond(getWorkflowState(workflowId).isTerminal)
    // Now that it's complete verify that we ended up in the state we're expecting
    getWorkflowState(workflowId) should equal (expectedState)
  }

  def getWorkflowState(workflowId: WorkflowId)(implicit ec: ExecutionContext): WorkflowState = {
    WorkflowState.fromString(getWorkflowMetadata(workflowId, Option(WorkflowMetadataKeys.Status)).getFields(WorkflowMetadataKeys.Status).headOption.map( _.asInstanceOf[JsString].value).getOrElse("Submitted"))
  }

  def getWorkflowOutputsFromMetadata(id: WorkflowId): Map[FullyQualifiedName, WdlValue] = {
    getWorkflowMetadata(id).fields.head._2.asInstanceOf[JsObject].getFields(WorkflowMetadataKeys.Outputs).toList match {
      case head::_ => head.asInstanceOf[JsObject].fields.map( x => (x._1, jsValueToWdlValue(x._2)))
      case _ => Map.empty
    }
  }

  private def jsValueToWdlValue(jsValue: JsValue): WdlValue = {
    jsValue match {
      case str: JsString => WdlString(str.value)
      case JsNumber(number) if number.scale == 0 => WdlInteger(number.intValue)
      case JsNumber(number) => WdlFloat(number.doubleValue)
      case JsBoolean(bool) => WdlBoolean(bool)
      case array: JsArray =>
        val valuesArray = array.elements.map(jsValueToWdlValue)
        if (valuesArray.isEmpty) WdlArray(WdlArrayType(WdlStringType), Seq.empty)
        else WdlArray(WdlArrayType(valuesArray.head.wdlType), valuesArray)
      case map: JsObject =>
        // TODO: currently assuming all keys are String. But that's not WDL-complete...
        val valuesMap: Map[WdlValue, WdlValue] = map.fields.map { case (fieldName, fieldValue) => (WdlString(fieldName), jsValueToWdlValue(fieldValue)) }
        if (valuesMap.isEmpty) WdlMap(WdlMapType(WdlStringType, WdlStringType), Map.empty)
        else WdlMap(WdlMapType(WdlStringType, valuesMap.head._2.wdlType), valuesMap)
    }
  }
}
