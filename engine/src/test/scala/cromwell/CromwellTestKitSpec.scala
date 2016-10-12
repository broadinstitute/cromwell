package cromwell

import java.nio.file.Paths
import java.util.UUID
import java.util.concurrent.atomic.AtomicInteger

import akka.actor.{Actor, ActorRef, ActorSystem, Props, Terminated}
import akka.pattern.ask
import akka.testkit._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestKitSpec._
import cromwell.backend._
import cromwell.core._
import cromwell.engine.backend.BackendConfigurationEntry
import cromwell.engine.workflow.WorkflowManagerActor.RetrieveNewWorkflows
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.{CacheLookupRequest, CacheResultMatchesForHashes}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.WorkflowSubmittedToStore
import cromwell.engine.workflow.workflowstore.{InMemoryWorkflowStore, WorkflowStoreActor}
import cromwell.jobstore.JobStoreActor.{JobStoreWriteSuccess, JobStoreWriterCommand}
import cromwell.server.{CromwellRootActor, CromwellSystem}
import cromwell.services.ServiceRegistryActor
import cromwell.services.metadata.MetadataQuery
import cromwell.services.metadata.MetadataService._
import cromwell.subworkflowstore.EmptySubWorkflowStoreActor
import cromwell.util.SampleWdl
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.metadata.MetadataBuilderActor
import org.scalactic.Equality
import org.scalatest.concurrent.{Eventually, ScalaFutures}
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, Matchers, OneInstancePerTest, WordSpecLike}
import spray.http.StatusCode
import spray.json._
import wdl4s.TaskCall
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}
import wdl4s.types._
import wdl4s.values._

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.matching.Regex

case class TestBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Set[TaskCall],
                                                serviceRegistryActor: ActorRef): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    throw new NotImplementedError("this is not implemented")
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
    NoFunctions
  }
}

case class OutputNotFoundException(outputFqn: String, actualOutputs: String) extends RuntimeException(s"Expected output $outputFqn was not found in: '$actualOutputs'")
case class LogNotFoundException(log: String) extends RuntimeException(s"Expected log $log was not found")

object CromwellTestKitSpec {
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
      |    # A dispatcher for actors performing blocking io operations
      |    # Prevents the whole system from being slowed down when waiting for responses from external resources for instance
      |    io-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |      # Using the forkjoin defaults, this can be tuned if we wish
      |    }
      |
      |    # A dispatcher for actors handling API operations
      |    # Keeps the API responsive regardless of the load of workflows being run
      |    api-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |
      |    # A dispatcher for engine actors
      |    # Because backends behaviour is unpredictable (potentially blocking, slow) the engine runs
      |    # on its own dispatcher to prevent backends from affecting its performance.
      |    engine-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |
      |    # A dispatcher used by supported backend actors
      |    backend-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |
      |    # Note that without further configuration, backend actors run on the default dispatcher
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
      |
      |services {}
    """.stripMargin

  val TimeoutDuration = 60 seconds

  private val testWorkflowManagerSystemCount = new AtomicInteger()

  class TestWorkflowManagerSystem extends CromwellSystem {
    override protected def systemName: String = "test-system-" + testWorkflowManagerSystemCount.incrementAndGet()
    override protected def newActorSystem() = ActorSystem(systemName, ConfigFactory.parseString(CromwellTestKitSpec.ConfigText))
    /**
      * Do NOT shut down the test actor system inside the normal flow.
      * The actor system will be externally shutdown outside the block.
      */
    // -Ywarn-value-discard
    override def shutdownActorSystem(): Future[Terminated] = { Future.successful(null) }

    def shutdownTestActorSystem() = super.shutdownActorSystem()
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

  /**
    * It turns out that tests using the metadata refresh actor running in parallel don't work so well if there are more
    * than one refresh actor. After many fits and starts the author decided that discretion was the better part of valor
    * as it is difficult to only singleton-ize the refresher or even the metadata service itself.
    *
    * Instead, doing things the "old way" for the tests, setting up a separate actor system and a singleton service
    * registry - this will be used by the CromwellRootActor within tests.
    *
    * :'(
    */
  private val ServiceRegistryActorSystem = akka.actor.ActorSystem("cromwell-service-registry-system")

  val ServiceRegistryActorInstance = {
    ServiceRegistryActorSystem.actorOf(ServiceRegistryActor.props(ConfigFactory.load()), "ServiceRegistryActor")
  }

  class TestCromwellRootActor(config: Config) extends CromwellRootActor {
    override val serverMode = true
    override lazy val serviceRegistryActor = ServiceRegistryActorInstance
    override lazy val workflowStore = new InMemoryWorkflowStore
    override val abortJobsOnTerminate = false
    def submitWorkflow(sources: WorkflowSourceFilesWithoutImports): WorkflowId = {
      val submitMessage = WorkflowStoreActor.SubmitWorkflow(sources)
      val result = Await.result(workflowStoreActor.ask(submitMessage)(TimeoutDuration), Duration.Inf).asInstanceOf[WorkflowSubmittedToStore].workflowId
      workflowManagerActor ! RetrieveNewWorkflows
      result
    }
  }
}

abstract class CromwellTestKitSpec(val twms: TestWorkflowManagerSystem = new CromwellTestKitSpec.TestWorkflowManagerSystem()) extends TestKit(twms.actorSystem)
  with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll with ScalaFutures with OneInstancePerTest with Eventually {

  override protected def afterAll() = { twms.shutdownTestActorSystem(); () }

  implicit val defaultPatience = PatienceConfig(timeout = Span(200, Seconds), interval = Span(1000, Millis))
  implicit val ec = system.dispatcher

  val dummyServiceRegistryActor = system.actorOf(Props.empty)
  val dummyLogCopyRouter = system.actorOf(Props.empty)

  // Allow to use shouldEqual between 2 WdlTypes while acknowledging for edge cases
  implicit val wdlTypeSoftEquality = new Equality[WdlType] {
    override def areEqual(a: WdlType, b: Any): Boolean = (a, b) match {
      case (WdlStringType | WdlFileType, WdlFileType | WdlStringType) => true
      case (arr1: WdlArrayType, arr2: WdlArrayType) => areEqual(arr1.memberType, arr2.memberType)
      case (map1: WdlMapType, map2: WdlMapType) => areEqual(map1.valueType, map2.valueType)
      case _ => a == b
    }
  }

  // Allow to use shouldEqual between 2 WdlValues while acknowledging for edge cases and checking for WdlType compatibilty
  implicit val wdlEquality = new Equality[WdlValue] {
    def fileEquality(f1: String, f2: String) = Paths.get(f1).getFileName == Paths.get(f2).getFileName

    override def areEqual(a: WdlValue, b: Any): Boolean = {
      val typeEquality = b match {
        case v: WdlValue => wdlTypeSoftEquality.areEqual(a.wdlType, v.wdlType)
        case _ => false
      }

      val valueEquality = (a, b) match {
        case (_: WdlFile, expectedFile: WdlFile) => fileEquality(a.valueString, expectedFile.valueString)
        case (_: WdlString, expectedFile: WdlFile) => fileEquality(a.valueString, expectedFile.valueString)
        case (array: WdlArray, expectedArray: WdlArray) =>
          (array.value.length == expectedArray.value.length) &&
            array.value.zip(expectedArray.value).map(Function.tupled(areEqual)).forall(identity)
        case (map: WdlMap, expectedMap: WdlMap) =>
          val mapped = map.value.map {
            case (k, v) => expectedMap.value.get(k).isDefined && areEqual(v, expectedMap.value(k))
          }

          (map.value.size == expectedMap.value.size) && mapped.forall(identity)
        case _ => a == b
      }

      typeEquality && valueEquality
    }
  }

  private def buildCromwellRootActor(config: Config) = {
    TestActorRef(new TestCromwellRootActor(config), name = "TestCromwellRootActor")
  }

  def runWdl(sampleWdl: SampleWdl,
             runtime: String = "",
             workflowOptions: String = "{}",
             terminalState: WorkflowState = WorkflowSucceeded,
             config: Config = DefaultConfig,
             patienceConfig: PatienceConfig = defaultPatience)(implicit ec: ExecutionContext): Map[FullyQualifiedName, WdlValue] = {
    val rootActor = buildCromwellRootActor(config)
    val sources = WorkflowSourceFilesWithoutImports(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, workflowOptions)
    val workflowId = rootActor.underlyingActor.submitWorkflow(sources)
    eventually { verifyWorkflowState(rootActor.underlyingActor.serviceRegistryActor, workflowId, terminalState) } (config = patienceConfig, pos = implicitly[org.scalactic.source.Position])
    val outcome = getWorkflowOutputsFromMetadata(workflowId, rootActor.underlyingActor.serviceRegistryActor)
    system.stop(rootActor)
    // And return the outcome:
    outcome
  }

  def runWdlAndAssertOutputs(sampleWdl: SampleWdl,
                             eventFilter: EventFilter,
                             expectedOutputs: Map[FullyQualifiedName, WdlValue],
                             runtime: String = "",
                             workflowOptions: String = "{}",
                             allowOtherOutputs: Boolean = true,
                             terminalState: WorkflowState = WorkflowSucceeded,
                             config: Config = DefaultConfig,
                             patienceConfig: PatienceConfig = defaultPatience)
                            (implicit ec: ExecutionContext): WorkflowId = {
    val rootActor = buildCromwellRootActor(config)
    val sources = sampleWdl.asWorkflowSources(runtime, workflowOptions)

    val workflowId = rootActor.underlyingActor.submitWorkflow(sources)
    eventually { verifyWorkflowState(rootActor.underlyingActor.serviceRegistryActor, workflowId, terminalState) } (config = patienceConfig, pos = implicitly[org.scalactic.source.Position])

    val outputs = getWorkflowOutputsFromMetadata(workflowId, rootActor.underlyingActor.serviceRegistryActor)
    val actualOutputNames = outputs.keys mkString ", "
    val expectedOutputNames = expectedOutputs.keys mkString " "

    expectedOutputs foreach { case (outputFqn, expectedValue) =>
      val actualValue = outputs.getOrElse(outputFqn, throw OutputNotFoundException(outputFqn, actualOutputNames))
      if (expectedValue != AnyValueIsFine) actualValue shouldEqual expectedValue
    }
    if (!allowOtherOutputs) {
      outputs foreach { case (actualFqn, actualValue) =>
        val expectedValue = expectedOutputs.getOrElse(actualFqn, throw new RuntimeException(s"Actual output $actualFqn was not wanted in '$expectedOutputNames'"))
        if (expectedValue != AnyValueIsFine) actualValue shouldEqual expectedValue
      }
    }

    system.stop(rootActor)
    workflowId
  }

  def getWorkflowMetadata(workflowId: WorkflowId, serviceRegistryActor: ActorRef, key: Option[String] = None)(implicit ec: ExecutionContext): JsObject = {
    // MetadataBuilderActor sends its response to context.parent, so we can't just use an ask to talk to it here
    val message = GetMetadataQueryAction(MetadataQuery(workflowId, None, key, None, None, expandSubWorkflows = false))
    val parentProbe = TestProbe()

    TestActorRef(MetadataBuilderActor.props(serviceRegistryActor), parentProbe.ref, s"MetadataActor-${UUID.randomUUID()}") ! message
    val metadata = parentProbe.expectMsgPF(TimeoutDuration) {
      // Because of type erasure the scala compiler can't check that the RequestComplete generic type will be (StatusCode, JsObject), which would generate a warning
      // As long as Metadata sends back a JsObject this is safe
      case response: RequestComplete[(StatusCode, JsObject)] @unchecked => response.response._2
    }

    system.stop(parentProbe.ref)
    metadata
  }

  /**
    * Verifies that a state is correct. // TODO: There must be a better way...?
    */
  private def verifyWorkflowState(serviceRegistryActor: ActorRef, workflowId: WorkflowId, expectedState: WorkflowState)(implicit ec: ExecutionContext): Unit = {
    def getWorkflowState(workflowId: WorkflowId, serviceRegistryActor: ActorRef)(implicit ec: ExecutionContext): WorkflowState = {
      val statusResponse = serviceRegistryActor.ask(GetStatus(workflowId))(TimeoutDuration).collect {
        case StatusLookupResponse(_, state) => state
        case f => throw new RuntimeException(s"Unexpected status response for $workflowId: $f")
      }
      Await.result(statusResponse, TimeoutDuration)
    }

    getWorkflowState(workflowId, serviceRegistryActor) should equal (expectedState)
    ()
  }

  private def getWorkflowOutputsFromMetadata(id: WorkflowId, serviceRegistryActor: ActorRef): Map[FullyQualifiedName, WdlValue] = {
    getWorkflowMetadata(id, serviceRegistryActor, None).getFields(WorkflowMetadataKeys.Outputs).toList match {
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

class AlwaysHappyJobStoreActor extends Actor {
  override def receive: Receive = {
    case x: JobStoreWriterCommand => sender ! JobStoreWriteSuccess(x)
  }
}

object AlwaysHappySubWorkflowStoreActor {
  def props: Props = Props(new EmptySubWorkflowStoreActor)
}

object AlwaysHappyJobStoreActor {
  def props: Props = Props(new AlwaysHappyJobStoreActor)
}

class EmptyCallCacheReadActor extends Actor {
  override def receive: Receive = {
    case CacheLookupRequest(CallCacheHashes(hashes)) => sender ! CacheResultMatchesForHashes(hashes, Set.empty)
  }
}

object EmptyCallCacheReadActor {
  def props: Props = Props(new EmptyCallCacheReadActor)
}
