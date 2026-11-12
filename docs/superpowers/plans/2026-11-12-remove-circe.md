# Remove circe (and jawn-parser) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move every circe usage onto spray-json so that no `io.circe:*` or `org.typelevel:jawn-*` artifact remains on any Cromwell classpath, eliminating CVE-2026-59990 by removal rather than by version override.

**Architecture:** spray-json is already Cromwell's primary JSON library (83 main files vs circe's 23). Each circe site is ported in place to spray-json `JsonWriter`/`JsonReader`/`RootJsonFormat` instances and `PrettyPrinter`/`CompactPrinter`. Libraries that only bridge to circe (`akka-http-circe`, `http4s-circe`, `sttp-circe`, `circe-config`, `circe-yaml`) are replaced by things already on the classpath (`akka-http-spray-json`, http4s `EntityDecoder.text`, ficus). YAML workflow-input support is dropped (user-approved breaking change, 2026-11-12); inputs become JSON-only.

**Tech Stack:** Scala 2.13.18, sbt, spray-json 1.3.6, akka-http-spray-json, http4s 0.21.31, ficus.

**Spec:** No separate spec doc. Requirements come from the 2026-11-12 session: "We want to remove vulnerable dependencies instead of updating them, when possible", "consolidate down to one [JSON parser]", "I'm OK with breaking YAML format input files".

## Global Constraints

- Done means: for `server`, `cromiam`, `centaur`, `cromwell-drs-localizer`, `womtool`, `wes2cromwell`, `export <p>/Runtime/dependencyClasspath` (and `Test/` for centaur) contains no `circe` and no `jawn` jar.
- No new third-party libraries. Use only what is already resolved (spray-json, akka-http-spray-json, ficus, http4s core).
- Smallest change that works: port each site in place; no shared "json utils" module, no drive-by cleanups.
- Error-response bodies from `WebServiceUtils.failRequest`/`errorRequest` stay **byte-identical** (customers display them raw — see comment in `WebServiceUtils.scala`).
- `/api/womtool/v1/describe` responses stay **semantically identical**, including explicit `"default": null` for inputs without defaults. Key order may change (spray `JsObject` is a `Map`).
- Internal `toString`/log JSON (execution store, value store, heartbeat config, token dispenser) need only remain valid, readable JSON; exact layout is not a contract.
- Every task ends with `sbt Test/compile` and `sbt assembly` passing (project rule), plus the task's own tests.
- After each task, `git grep -n "io.circe\|akkahttpcirce\|org.http4s.circe" -- <task's module paths>` returns nothing.

## Review Focus

1. **YAML-only syntax in inputs** (`# comments`, unquoted keys, single quotes, trailing commas): should be rejected with HTTP 400 and the "not valid JSON" message, not a 500. → Task 2 `rejectYamlOnlySyntax`.
2. **Deeply nested JSON inputs** (10,000 nested `[`), the CVE's attack class: should be a recoverable 400, not `StackOverflowError`. → Task 2 `rejectDeepNesting`.
3. **Blank / whitespace-only inputs string**: should behave exactly as on `develop` today (record the current result in a characterization test before changing code). → Task 2 `blankInputsUnchanged`.
4. **Big or precise numbers in inputs** (`12345678901234567890`, `0.1000`, `1e10`): should survive into the stored inputs JSON without precision loss. → Task 2 `numbersRoundTrip`.
5. **Error messages with quotes, newlines, non-ASCII, and multi-element `errors` arrays**: should print byte-identically to today. → Task 4 `errorResponseGolden`.

---

### Task 1: Drop circe-only transitive dependencies (gating task)

`workbench-google` (in `cloudSupportDependencies`) pulls `workbench-util2` → `circe-fs2` → `circe-jawn` → `jawn-parser` into nearly every project, and no Cromwell source imports a workbench class outside CromIAM (which uses the separate `workbench-model`/`workbench-util`). `sttp-circe` in `drsLocalizerDependencies` has no sttp usage in `cromwell-drs-localizer`. If `workbench-google` turns out to be needed, **stop and escalate**: jawn cannot be removed without it gone, and the remaining tasks lose their purpose.

**Files:**
- Modify: `project/Dependencies.scala:298-299` (remove `workbench-google` entry), `project/Dependencies.scala:494` (remove `"com.softwaremill.sttp" %% "circe"`), `project/Dependencies.scala:102` (remove `workbenchGoogleV` if now unused)

**Interfaces:** Produces: no `circe-fs2`, `circe-jawn`, or `workbench-util2` on any classpath.

- [ ] **Step 1: Snapshot classpaths** — `sbt -batch "export server/Runtime/dependencyClasspath" "export cromwell-drs-localizer/Runtime/dependencyClasspath" "export cromiam/Runtime/dependencyClasspath"`, save one jar name per line, sorted.
- [ ] **Step 2: Remove the two dependencies.**
- [ ] **Step 3: Compile and assemble** — `sbt Test/compile` then `sbt assembly`. Expected: both succeed.
- [ ] **Step 4: Bytecode scan** — for each assembly jar, `unzip -l` and confirm no `org/broadinstitute/dsde/workbench/google/` or `workbench/util2/` classes remain, then `javap`-scan (or `grep -c` the extracted class files) for references to `org/broadinstitute/dsde/workbench/google` from Cromwell classes. Expected: zero references.
- [ ] **Step 5: Diff classpaths** — re-export and diff against Step 1. Expected: only `workbench-google`, `workbench-util2`, `circe-fs2`, `circe-jawn`, `sttp circe` and their exclusive transitives disappear; no other library changes version. Record the removed list in the commit message.
- [ ] **Step 6: Run** `sbt cloudSupport/test cromwell-drs-localizer/test`. Expected: PASS.
- [ ] **Step 7: Commit** `"[CTM-xxx] Remove unused workbench-google and sttp-circe dependencies"`.

---

### Task 2: JSON-only workflow inputs (drop YAML)

**Files:**
- Modify: `engine/src/main/scala/cromwell/webservice/PartialWorkflowSources.scala:200-246`
- Delete: `wom/src/main/scala/wom/util/YamlUtils.scala`, `wom/src/test/scala/wom/util/YamlUtilsSpec.scala`, the `yaml { ... }` block in `wom/src/main/resources/reference.conf`
- Modify: `docs/Configuring.md` (remove the `### YAML` section, ~lines 785-816), `CHANGELOG.md` (new entry under the unreleased version: "Workflow inputs must be JSON; YAML inputs are no longer accepted. The `yaml.max-nodes` and `yaml.max-depth` settings are removed.")
- Test: `engine/src/test/scala/cromwell/webservice/PartialWorkflowSourcesSpec.scala`

**Interfaces:**
- Produces: `private def arrayTypeElementValidation[A](data: String, interpretElementFunction: JsValue => ErrorOr[A]): ErrorOr[Vector[A]]` using `spray.json.JsonParser(data)`; `workflowInputsValidation(data: String): ErrorOr[Vector[WorkflowJson]]` (signature unchanged) prints each element with `compactPrint`; `workflowIdsValidation` signature unchanged, matching `JsString`.
- Error message becomes: `s"Input file is not valid JSON. Inputs data: '$data'. Error: ${ExceptionUtils.getMessage(error)}."` (same shape, "a valid yaml or json" → "valid JSON").

- [ ] **Step 1: Characterize blank input on develop (before any edit)** — add `blankInputsUnchanged`: call `workflowInputsValidation("  ")` and `workflowInputsValidation("")`, run it, and hard-code the observed results as the assertions. Commit this test alone so the record precedes the change.
- [ ] **Step 2: Write the failing tests**

```scala
it should "rejectYamlOnlySyntax" in {
  List("a: 1", "{'a': 1}", "{a: 1}", "{\"a\": 1,}", "# c\n{\"a\": 1}") foreach { yamlish =>
    val result = PartialWorkflowSources.workflowInputsValidation(yamlish)
    result.isInvalid shouldBe true
    result.swap.toOption.get.head should startWith("Input file is not valid JSON.")
  }
}

it should "rejectDeepNesting" in {
  val deep = "[" * 10000 + "]" * 10000
  val result = PartialWorkflowSources.workflowInputsValidation(deep) // must not throw StackOverflowError
  result.isInvalid shouldBe true
}

it should "numbersRoundTrip" in {
  val input = """{"w.big":12345678901234567890,"w.prec":0.1000,"w.exp":1e10}"""
  val Valid(Vector(out)) = PartialWorkflowSources.workflowInputsValidation(input)
  out.parseJson shouldBe input.parseJson
  out should include("12345678901234567890")
  out should include("0.1000")
}
```

Keep the existing JSON-object and JSON-array cases in the spec (lines ~60-120); they must still pass unchanged.

- [ ] **Step 3: Run** `sbt "engine/testOnly cromwell.webservice.PartialWorkflowSourcesSpec"`. Expected: `rejectYamlOnlySyntax` FAILS (YAML still accepted); others may pass.
- [ ] **Step 4: Implement** — replace `YamlUtils.parse` with `Try(JsonParser(data))`; `JsArray(elements)` → traverse elements, otherwise wrap single value. If `rejectDeepNesting` throws `StackOverflowError`, pass `JsonParserSettings.default.withMaxDepth(1000)` explicitly (spray-json 1.3.6 supports it); do not catch `StackOverflowError`.
- [ ] **Step 5: Delete** `YamlUtils`, `YamlUtilsSpec`, the `yaml` config block; update docs and CHANGELOG.
- [ ] **Step 6: Run** the spec, then `sbt Test/compile` and `sbt assembly`. Expected: PASS. `git grep -n YamlUtils` returns nothing.
- [ ] **Step 7: Commit** `"[CTM-xxx] Accept only JSON workflow inputs; drop YAML parsing"`.

---

### Task 3: Womtool describe models to spray-json

**Files:**
- Modify: `services/src/main/scala/cromwell/services/womtool/models/{InputDescription,OutputDescription,WomTypeJsonSupport,WomExpressionJsonSupport,MetaValueElementJsonSupport,WorkflowDescription}.scala`
- Test: `services/src/test/scala/cromwell/services/womtool/DescriberSpec.scala` (parse fixtures with `.parseJson`, compare `description.toJson shouldBe expected`)

**Interfaces:**
- Produces (encoder implicit names are kept, except `WorkflowDescription`'s encoder/decoder pair, which becomes one format):
  - `WomTypeJsonSupport.womTypeEncoder: JsonWriter[WomType]` — same object shapes as today.
  - `WomExpressionJsonSupport.womExpressionEncoder: JsonWriter[WomExpression]` — `JsString(sourceString)`.
  - `MetaValueElementJsonSupport.metaValueElementEncoder: JsonWriter[MetaValueElement]` — floats via `JsNumber(Double)`.
  - `InputDescription.inputDescriptionWriter: JsonWriter[InputDescription]` — hand-written, always emits all five fields, `default` as `JsNull` when `None`.
  - `OutputDescription.outputDescriptionWriter: JsonWriter[OutputDescription]`.
  - `WorkflowDescription.workflowDescriptionFormat: RootJsonFormat[WorkflowDescription]` — `write` emits all 12 fields; `read` keeps today's partial decode of `valid`, `errors`, `validWorkflow` only (needed by `WomtoolRouteSupportSpec`). Replaces both `workflowDescriptionEncoder` and `workflowDescriptionDecoder`.

- [ ] **Step 1: Port `DescriberSpec` to spray** (fixture parsing + `check` using `toJson`) while models are still circe, so it fails to compile. That is the failing test; the 29 `description.json` fixtures are the golden data and must not be edited.
- [ ] **Step 2: Run** `sbt "services/testOnly cromwell.services.womtool.DescriberSpec"`. Expected: compile FAIL (no `JsonWriter[WorkflowDescription]`).
- [ ] **Step 3: Implement the writers/format listed above.**
- [ ] **Step 4: Run** the spec. Expected: PASS (30 tests). Any fixture mismatch is a bug in the port, never a fixture edit.
- [ ] **Step 5: Commit** `"[CTM-xxx] Womtool describe models use spray-json"`.

---

### Task 4: Engine web layer (describe route, error responses)

**Files:**
- Modify: `engine/src/main/scala/cromwell/webservice/routes/WomtoolRouteSupport.scala:88-92`, `engine/src/main/scala/cromwell/webservice/WebServiceUtils.scala`
- Test: create `engine/src/test/scala/cromwell/webservice/WebServiceUtilsSpec.scala`; modify `engine/src/test/scala/cromwell/webservice/routes/WomtoolRouteSupportSpec.scala` (replace `ErrorAccumulatingCirceSupport._` with `SprayJsonSupport._`)

**Interfaces:**
- Consumes: `WorkflowDescription.workflowDescriptionFormat` (Task 3); existing `WorkflowJsonSupport.errorResponse: RootJsonFormat[FailureResponse]`.
- Produces: `WebServiceUtils.prettyPrint(failureResponse: FailureResponse): String` (still private) built on a private `object CirceLayoutPrettyPrinter extends PrettyPrinter` overriding `printArray` to emit one element per line, indented, matching circe `Printer.spaces2.copy(dropNullValues = true, colonLeft = "")`. Describe route does `complete(response.description)` via `SprayJsonSupport`.

- [ ] **Step 1: Write the golden test on develop's circe code first**

```scala
it should "errorResponseGolden" in {
  val ex = AggregatedMessageException("Bad \"thing\"\nhappened: é", List("first", "second \"quoted\""))
  Get() ~> ex.failRequest(StatusCodes.BadRequest) ~> check {
    responseAs[String] shouldBe
      """{
        |  "status": "fail",
        |  "message": "Bad \"thing\"\nhappened: é",
        |  "errors": [
        |    "first",
        |    "second \"quoted\""
        |  ]
        |}""".stripMargin
  }
  Get() ~> new RuntimeException("solo").errorRequest(StatusCodes.InternalServerError) ~> check {
    responseAs[String] shouldBe "{\n  \"status\": \"error\",\n  \"message\": \"solo\"\n}"
  }
}
```

Run it against the unmodified circe implementation and adjust literals only if circe's real output differs from the above (the test records today's bytes). Commit the test alone.
- [ ] **Step 2: Port** `WebServiceUtils` (drop `io.circe._`, `io.circe.generic.auto._`) and `WomtoolRouteSupport` as above; port `WomtoolRouteSupportSpec`.
- [ ] **Step 3: Run** `sbt "engine/testOnly cromwell.webservice.WebServiceUtilsSpec cromwell.webservice.routes.WomtoolRouteSupportSpec cromwell.webservice.routes.CromwellApiServiceSpec"`. Expected: PASS, golden unchanged.
- [ ] **Step 4: Commit** `"[CTM-xxx] Engine web layer uses spray-json; drop akka-http-circe"`.

---

### Task 5: Engine/core internal JSON (logging and toString)

**Files:**
- Modify: `core/src/main/scala/cromwell/core/JobToken.scala`, `engine/src/main/scala/cromwell/engine/workflow/tokens/{JobTokenDispenserActor,TokenQueue,UnhoggableTokenPool}.scala`, `engine/src/main/scala/cromwell/engine/workflow/lifecycle/execution/stores/{ExecutionStore,ValueStore}.scala`, `engine/src/main/scala/cromwell/engine/workflow/workflowstore/WorkflowHeartbeatConfig.scala`
- Test: `engine/src/test/scala/cromwell/engine/workflow/workflowstore/WorkflowHeartbeatConfigSpec.scala` (exists or create)

**Interfaces:**
- Produces: `JobTokenType` companion `implicit val jobTokenTypeFormat: RootJsonFormat[JobTokenType] = jsonFormat3(...)` (in `core`); `jsonFormatN` formats for `HogGroupState`, `TokenPoolState`, `TokenQueueSize`, `TokenQueueState`, `TokenTypeState`, `TokenDispenserState` in their companions; `WorkflowHeartbeatConfig` keeps a custom `JsonWriter[FiniteDuration]` writing `JsString(d.toString)`. All `toString`s use `.toJson.prettyPrint`.

- [ ] **Step 1: Characterization test** — `heartbeatToStringIsJson`: `WorkflowHeartbeatConfig(ConfigFactory.load()).toString.parseJson.asJsObject.fields("ttl") shouldBe JsString("10 minutes")` (use the reference.conf default ttl, read it from config in the test rather than hard-coding if it differs).
- [ ] **Step 2: Run** it on the circe code. Expected: PASS. This is a characterization guard (the output is a plain string either way); it must still pass after the port.
- [ ] **Step 3: Port** all listed files.
- [ ] **Step 4: Run** `sbt "engine/testOnly *WorkflowHeartbeatConfig* *JobTokenDispenser* *TokenQueue* *UnhoggableTokenPool* *ExecutionStore* *ValueStore*" core/test`. Expected: PASS.
- [ ] **Step 5: Commit** `"[CTM-xxx] Token, store, and heartbeat JSON use spray-json"`.

---

### Task 6: DRS resolver

**Files:**
- Modify: `cloud-nio/cloud-nio-impl-drs/src/main/scala/cloud/nio/impl/drs/DrsPathResolver.scala:11-14,63,118,218-300`
- Test: `cloud-nio/cloud-nio-impl-drs/src/test/scala/cloud/nio/impl/drs/DrsPathResolverSpec.scala`

**Interfaces:**
- Produces: `SADataObject(data: JsValue)` (was circe `Json`; only tests construct it). In the resolver's companion: `JsonWriter[DrsResolverField.Value]` and `JsonWriter[DrsCloudPlatform.Value]` writing `JsString(v.toString)`; `RootJsonWriter[DrsResolverRequest]`; `RootJsonReader`s for `DrsResolverResponse`, `AccessUrl`, `SADataObject`, `DrsResolverFailureResponse`, `DrsResolverFailureResponsePayload`. Readers must ignore unknown fields and treat missing **or null** optional fields as `None` (spray `jsonFormatN` does both).
- Request body: `makeDrsResolverRequest(...).toJson.compactPrint`. Response: `Try(responseContent.parseJson.convertTo[DrsResolverResponse])` lifted to `IO`.

- [ ] **Step 1: Failing test** — `decodeNullsAndUnknownFields`: `"""{"size":null,"bogus":1,"gsUri":"gs://b/o"}"""` decodes to `DrsResolverResponse(size = None, gsUri = Some("gs://b/o"))`. Port `mockGSA` to `JsObject("key" -> JsString("value"))`.
- [ ] **Step 2: Run** `sbt "cloud-nio-impl-drs/testOnly cloud.nio.impl.drs.DrsPathResolverSpec"`. Expected: compile FAIL.
- [ ] **Step 3: Port** the resolver.
- [ ] **Step 4: Run** the spec plus `sbt drsFileSystem/test cromwell-drs-localizer/test`. Expected: PASS.
- [ ] **Step 5: Commit** `"[CTM-xxx] DRS resolver uses spray-json"`.

---

### Task 7: Docker registry hashing

**Files:**
- Modify: `dockerHashing/src/main/scala/cromwell/docker/registryv2/DockerRegistryV2Abstract.scala:10-14,59-80`; formats in `dockerHashing/src/main/scala/cromwell/docker/registryv2/{DockerManifest,DockerAccessToken}.scala` companions
- Modify: `project/Dependencies.scala:131` (remove `http4s-circe`), `:716` (`val http4sCirce`, unused in build.sbt)
- Test: `dockerHashing/src/test/scala/cromwell/docker/registryv2/DockerRegistryV2AbstractSpec.scala` (exists or create)

**Interfaces:**
- Produces: `private def sprayEntityDecoder[A: JsonReader](mediaRange: MediaRange): EntityDecoder[IO, A]` built from `EntityDecoder.decodeBy(mediaRange)` + `EntityDecoder.decodeText`, mapping parse/convert failures to `MalformedMessageBodyFailure`. Token decoder uses it with `MediaType.application.json`. Formats: `jsonFormat1` for `DockerLayer`, `DockerManifest`, `DockerManifestReference`, `DockerManifestList`, `DockerAccessToken` (Array fields need `ClassTag`, already available).

- [ ] **Step 1: Failing test** — `decodeVendorManifest`: a `Response[IO]` with `Content-Type: application/vnd.docker.distribution.manifest.v2+json` and body `{"schemaVersion":2,"layers":[{"size":10},{"size":32}]}` decodes via `entityManifestDecoder` to a manifest with `compressedSize == 42`; a body of `not json` yields `Left(MalformedMessageBodyFailure)`.
- [ ] **Step 2: Run** `sbt "dockerHashing/testOnly *DockerRegistryV2Abstract*"`. Expected: FAIL to compile or FAIL on the malformed case.
- [ ] **Step 3: Port** and remove `http4s-circe`.
- [ ] **Step 4: Run** `sbt dockerHashing/test`. Expected: PASS (33 + new).
- [ ] **Step 5: Commit** `"[CTM-xxx] Docker hashing decodes manifests with spray-json; drop http4s-circe"`.

---

### Task 8: GCP Batch backend

**Files:**
- Modify: `supportedBackends/google/batch/src/main/scala/cromwell/backend/google/batch/actors/GcpBatchInitializationActor.scala:3-5,170-171`; `supportedBackends/google/batch/src/main/scala/cromwell/backend/google/batch/models/GcpBatchConfigurationAttributes.scala:399-417`
- Test: `supportedBackends/google/batch/src/test/scala/cromwell/backend/google/batch/models/GcpBatchConfigurationAttributesSpec.scala`

**Interfaces:**
- Produces: `ProjectLabels` decoded with `jsonFormat1` (`Try(httpResponse.parseAsString().parseJson.convertTo[ProjectLabels])`). Reference-disk manifests read with ficus `config.as[ManifestFile]` via `net.ceedubs.ficus.readers.ArbitraryTypeReader._` (field names map to keys as-is, same as circe-config), wrapped in `Try` per config; failure message keeps the prefix `s"Reference disks misconfigured for backend $backendName, could not parse as List[ManifestFile]"`. Delete the comment about the `as` collision.

- [ ] **Step 1: Failing test** — `referenceDiskManifestMissingField`: a reference-disk config entry without `diskSizeGb` produces an invalid result whose message starts with the prefix above; a valid entry with `crc32c = 4294967295` (exceeds Int) parses with that `Long` value. Confirm existing reference-disk tests in the spec still describe today's accepted configs.
- [ ] **Step 2: Run** `sbt "googleBatch/testOnly *GcpBatchConfigurationAttributesSpec*"` (check exact project ID in `build.sbt`). Expected: the new test FAILS or passes on circe (keep it either way as a guard).
- [ ] **Step 3: Port** both files.
- [ ] **Step 4: Run** the batch backend's full test suite. Expected: PASS.
- [ ] **Step 5: Commit** `"[CTM-xxx] GCP Batch backend uses spray-json and ficus instead of circe"`.

---

### Task 9: Centaur

**Files:**
- Modify: `centaur/src/main/scala/centaur/CromwellTracker.scala`, `centaur/src/main/scala/centaur/test/metadata/CallAttemptFailure.scala`, `centaur/src/main/scala/centaur/test/Test.scala:35,336-339`
- Test: `centaur/src/test/scala/centaur/test/metadata/CallAttemptFailureSpec.scala`

**Interfaces:**
- Produces: `CallAttemptFailure.buildFailures(json: String): IO[Vector[CallAttemptFailure]]` (signature unchanged). Keep current leniency exactly: a malformed top level yields `Vector.empty`; a malformed or incomplete attempt (missing `shardIndex`, `attempt`, or `failures[0].message`) is skipped, not an error. `CromwellTracker` reads `workflowProcessingEvents` from the parsed `JsObject` (replacing circe-optics `root....each`) with a `jsonFormat` for `WorkflowProcessingEvent`; `toString` uses `prettyPrint`. `Test.scala` uses `parseJson`/`asJsObject.fields.get("failures")`/`prettyPrint`.

- [ ] **Step 1: Failing test** — add `skipIncompleteAttempts` to `CallAttemptFailureSpec`: metadata with one call whose attempts are `[{"shardIndex":0,"attempt":1,"failures":[{"message":"m"}]}, {"attempt":1}]` yields exactly one failure; `"not json"` yields `Vector.empty`.
- [ ] **Step 2: Run** `sbt "centaur/testOnly centaur.test.metadata.CallAttemptFailureSpec"`. Expected: new test passes on circe (characterization). Keep it.
- [ ] **Step 3: Port** the three files.
- [ ] **Step 4: Run** `sbt centaur/test`. Expected: PASS (existing `parse metadata` test still returns 10).
- [ ] **Step 5: Commit** `"[CTM-xxx] Centaur uses spray-json"`.

---

### Task 10: Womtool outputs, build cleanup, and final verification

**Files:**
- Modify: `womtool/src/main/scala/womtool/outputs/Outputs.scala` (spray `JsObject`/`JsString`, `prettyPrint`)
- Modify: `project/Dependencies.scala`: delete `circeV`, `circeYamlV`, `circeConfigV`, `akkaHttpCirceIntegrationV`, `circeYamlDependency`, `circeDependencies` and every `++ circeDependencies` (lines ~113-119, 169, 390, 416, 423, 434, 440, 498), `akkaHttpCirceIntegrationDependency` and its use at 440, and the `jawnOverrides` list plus its entry in `cromwellDependencyOverrides` (added 2026-11-12; jawn no longer resolves).
- Test: existing `womtool/src/test` outputs tests.

- [ ] **Step 1: Port `Outputs.scala`**, then run `sbt "womtool/testOnly *Outputs*"` (or `womtool/test`). Expected: PASS.
- [ ] **Step 2: Remove circe from the build** as listed.
- [ ] **Step 3: Verify no source references** — `git grep -n "io\.circe\|akkahttpcirce\|org\.http4s\.circe\|jawn" -- '*.scala'`. Expected: no output.
- [ ] **Step 4: Verify classpaths** — export `Runtime/dependencyClasspath` for `server`, `cromiam`, `womtool`, `cromwell-drs-localizer`, `wes2cromwell` and `Test/dependencyClasspath` for `centaur`; `grep -i "circe\|jawn"`. Expected: no output. Diff against Task 1's snapshot: only circe/jawn/akka-http-circe/http4s-circe/http4s-jawn/jawn-fs2/circe-yaml(+its snakeyaml only if nothing else needs it) disappear.
- [ ] **Step 5: Verify assemblies** — `sbt assembly`, then `unzip -l` each assembly and `grep -c "io/circe\|org/typelevel/jawn"`. Expected: 0.
- [ ] **Step 6: Full checks** — `sbt Test/compile`, `sbt test`. Expected: PASS.
- [ ] **Step 7: Commit** `"[CTM-xxx] Remove circe and jawn from the build"`.
