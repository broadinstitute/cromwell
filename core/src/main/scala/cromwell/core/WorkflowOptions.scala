package cromwell.core

import cats.syntax.validated._
import com.typesafe.config.ConfigFactory
import common.util.TryUtil
import common.validation.ErrorOr.ErrorOr
import spray.json._

import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

/**
 * WorkflowOptions is constructed with a String -> String map (usually as a JsObject).
 *
 * During construction, any keys in the map that are in the workflow-options.encrypted-fields list
 * will be replaced with a JSON object which represents the encrypted value. For example,
 * if the initial map is:
 *
 * {
 *   "encrypt_this": "foo",
 *   "dont_encrypt": "bar"
 * }
 *
 * Then upon construction, if workflow-options.encrypted-fields=["encrypt_this"], then the resulting
 * JSON object is:
 *
 * {
 *   "encrypt_this": {
 *       "iv": "... base64 encoded initialization vector ...",
 *       "ciphertext": "... base64 encoded cipher text ..."
 *   },
 *   "dont_encrypt": "bar"
 * }
 *
 * When the client calls workflowOptions.get("encrypt_this"), the get method will recognize that the
 * value is encrypted and decrypt it before sending back the unencrypted value.
 *
 * The reason why there's a .fromJsonObject(obj: JsObject), is that we need to be able to construct
 * workflow options from JSON that looks like the encrypted example above, like if restart the server
 * and we need to reload what's in the database back into a WorkflowOptions object.
 *
 * There's logic in here that says that WorkflowOptions can only be constructed from JsObject if the
 * values are either JsString or a JsObject that looks like {"iv": "...", "ciphertext": "..."}
 */

object WorkflowOptions {

  sealed abstract class WorkflowOption(val name: String)
  // Caching
  case object WriteToCache extends WorkflowOption("write_to_cache")
  case object ReadFromCache extends WorkflowOption("read_from_cache")
  case object CallCacheHitPathPrefixes extends WorkflowOption("call_cache_hit_path_prefixes")

  // Copying
  case object FinalWorkflowLogDir extends WorkflowOption("final_workflow_log_dir")
  case object FinalCallLogsDir extends WorkflowOption("final_call_logs_dir")
  case object FinalWorkflowOutputsDir extends WorkflowOption("final_workflow_outputs_dir")
  case object UseRelativeOutputPaths extends WorkflowOption(name = "use_relative_output_paths")
  case object FinalWorkflowOutputsDirMetadata extends WorkflowOption("final_workflow_outputs_dir_metadata") {
    // Default to Source because that was originally the only behavior
    def fromString(s: Option[String]): FinalWorkflowOutputsDirMetadata =
      s match {
        case Some("source") => Source
        case Some("destination") => Destination
        case _ => Source
      }
  }

  sealed trait FinalWorkflowOutputsDirMetadata
  // Metadata points to original location (copy source)
  case object Source extends FinalWorkflowOutputsDirMetadata
  // Metadata points to new location (copy destination)
  case object Destination extends FinalWorkflowOutputsDirMetadata

  // Misc.
  case object DefaultRuntimeOptions extends WorkflowOption("default_runtime_attributes")
  case object WorkflowFailureMode extends WorkflowOption("workflow_failure_mode")
  case object UseReferenceDisks extends WorkflowOption("use_reference_disks")
  case object MemoryRetryMultiplier extends WorkflowOption("memory_retry_multiplier")
  case object WorkflowCallbackUri extends WorkflowOption("workflow_callback_uri")

  private lazy val WorkflowOptionsConf = ConfigFactory.load.getConfig("workflow-options")
  private lazy val EncryptedFields: Seq[String] = WorkflowOptionsConf.getStringList("encrypted-fields").asScala.toList
  private lazy val EncryptionKey: String = WorkflowOptionsConf.getString("base64-encryption-key")
  private lazy val defaultRuntimeOptionKey: String = DefaultRuntimeOptions.name
  private lazy val validObjectKeys: Set[String] = Set(DefaultRuntimeOptions.name, "google_labels")

  def encryptField(value: JsString): Try[JsObject] =
    Aes256Cbc.encrypt(value.value.getBytes("utf-8"), SecretKey(EncryptionKey)) match {
      case Success(encryptedValue) =>
        Success(
          JsObject(
            Map(
              "iv" -> JsString(encryptedValue.base64Iv),
              "ciphertext" -> JsString(encryptedValue.base64CipherText)
            )
          )
        )
      case Failure(ex) => Failure(ex)
    }

  def decryptField(obj: JsObject): Try[String] =
    (obj.fields.get("iv"), obj.fields.get("ciphertext")) match {
      case (Some(iv: JsString), Some(ciphertext: JsString)) =>
        Aes256Cbc
          .decrypt(EncryptedBytes(ciphertext.value, iv.value), SecretKey(WorkflowOptions.EncryptionKey))
          .map(new String(_, "utf-8"))
      case _ => Failure(new RuntimeException(s"JsObject must have 'iv' and 'ciphertext' fields to decrypt: $obj"))
    }

  def isEncryptedField(jsValue: JsValue): Boolean = jsValue match {
    case obj: JsObject if obj.fields.keys.exists(_ == "iv") && obj.fields.keys.exists(_ == "ciphertext") => true
    case _ => false
  }

  def fromJsonObject(jsObject: JsObject): Try[WorkflowOptions] = {
    val encrypted: Map[String, Try[JsValue]] = jsObject.fields map {
      case (k, v: JsString) if EncryptedFields.contains(k) => k -> encryptField(v)
      case (k, v: JsString) => k -> Success(v)
      case (k, v: JsBoolean) => k -> Success(v)
      case (k, v: JsObject) if validObjectKeys.contains(k) => k -> Success(v)
      case (k, v: JsNumber) => k -> Success(v)
      case (k, v) if isEncryptedField(v) => k -> Success(v)
      case (k, v: JsArray) => k -> Success(v)
      case (k, v) =>
        k -> Failure(new UnsupportedOperationException(s"Unsupported key/value pair in WorkflowOptions: $k -> $v"))
    }

    encrypted.values collect { case f: Failure[_] => f } match {
      case s if s.nonEmpty => Failure(new RuntimeException(s.map(_.exception.getMessage).mkString("\n")))
      case _ => Success(WorkflowOptions(new JsObject(encrypted map { case (k, v) => k -> v.get })))
    }
  }

  def fromJsonString(json: String): Try[WorkflowOptions] = Try(json.parseJson) match {
    case Success(obj: JsObject) => fromJsonObject(obj)
    case Failure(ex) => Failure(ex)
    case Success(x) => Failure(new UnsupportedOperationException(s"Expecting JSON object, got $x"))
  }

  def fromMap(m: Map[String, String]) = fromJsonObject(JsObject(m map { case (k, v) => k -> JsString(v) }))

  private def getAsJson(key: String, jsObject: JsObject) = jsObject.fields.get(key) match {
    case Some(jsStr: JsString) => Success(jsStr)
    case Some(jsNum: JsNumber) => Success(jsNum)
    case Some(jsArray: JsArray) => Success(jsArray)
    case Some(jsBool: JsBoolean) => Success(jsBool)
    case Some(jsObj: JsObject) if isEncryptedField(jsObj) => decryptField(jsObj) map JsString.apply
    case Some(jsObj: JsObject) => Success(jsObj)
    case Some(jsVal: JsValue) => Failure(new IllegalArgumentException(s"Unsupported value as JsValue: $jsVal"))
    case None => Failure(OptionNotFoundException(s"Field not found: $key"))
  }

  private def get(key: String, jsObject: JsObject) = jsObject.fields.get(key) match {
    case Some(jsStr: JsString) => Success(jsStr.value)
    case Some(jsNum: JsNumber) => Success(jsNum.value.toString)
    case Some(jsBool: JsBoolean) => Success(jsBool.value.toString)
    case Some(jsObj: JsObject) if isEncryptedField(jsObj) => decryptField(jsObj)
    case Some(jsVal: JsValue) => Failure(new IllegalArgumentException(s"Unsupported value as JsValue: $jsVal"))
    case None => Failure(OptionNotFoundException(s"Field not found: $key"))
  }

  val empty = WorkflowOptions.fromMap(Map.empty).get
}

case class WorkflowOptions(jsObject: JsObject) {
  import WorkflowOptions._

  def toMap = jsObject.fields

  def get(key: String): Try[String] = WorkflowOptions.get(key, jsObject)
  def get(option: WorkflowOption): Try[String] = get(option.name)

  def getBoolean(key: String): Try[Boolean] = jsObject.fields.get(key) match {
    case Some(jsBool: JsBoolean) => Success(jsBool.value)
    case Some(jsVal: JsValue) => Failure(new IllegalArgumentException(s"Unsupported JsValue as JsBoolean: $jsVal"))
    case None => Failure(OptionNotFoundException(s"Field not found: $key"))
  }

  def getVectorOfStrings(key: String): ErrorOr[Option[Vector[String]]] = jsObject.fields.get(key) match {
    case Some(jsArr: JsArray) => Option(jsArr.elements collect { case e: JsString => e.value }).validNel
    case Some(jsVal: JsValue) => s"Unsupported JsValue as JsArray: $jsVal".invalidNel
    case _ => None.validNel
  }

  def getDefaultRuntimeOption(key: String): Try[JsValue] = defaultRuntimeOptions map { attributes =>
    attributes.getOrElse(key, throw OptionNotFoundException(s"Field not found $key"))
  }

  lazy val defaultRuntimeOptions = jsObject.fields.get(defaultRuntimeOptionKey) match {
    case Some(jsObj: JsObject) =>
      TryUtil.sequenceMap(jsObj.fields map { case (k, _) => k -> WorkflowOptions.getAsJson(k, jsObj) })
    case Some(jsVal) =>
      Failure(
        new IllegalArgumentException(
          s"Unsupported JsValue for $defaultRuntimeOptionKey: $jsVal. Expected a JSON object."
        )
      )
    case None => Failure(OptionNotFoundException(s"Cannot find definition for default runtime attributes"))
  }

  def getOrElse[B >: String](key: String, default: => B): B = get(key) match {
    case Success(value) => value
    case _ => default
  }

  override def toString: String = asPrettyJson
  def asPrettyJson: String = jsObject.prettyPrint

  /**
    * Returns a JSON representation of these workflow options where the encrypted values
    * have been replaced by the string "cleared".
    *
    * Used to protect encrypted values from being stored in metadata
    */
  def clearEncryptedValues: WorkflowOptions = {

    def revoke(o: JsObject): JsObject = {
      val newFields = o.fields map {
        case (k, v: JsObject) if isEncryptedField(v) => k -> JsString("cleared")
        case (k, v: JsObject) => k -> revoke(v)
        case (k, v) => k -> v
      }
      JsObject(newFields)
    }

    WorkflowOptions(revoke(jsObject))
  }
}

case class OptionNotFoundException(message: String) extends Exception(message)
