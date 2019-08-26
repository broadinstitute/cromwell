package cromwell

import cromwell.CromwellApp.{Run, Server, Submit}
import cromwell.CromwellCommandLineSpec.WdlAndInputs
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{FileClobber, FilePassingWorkflow, ThreeStep}
import org.scalatest.{BeforeAndAfter, FlatSpec, Matchers}

import scala.util.Try

class CromwellCommandLineSpec extends FlatSpec with Matchers with BeforeAndAfter {

  var parser: scopt.OptionParser[CommandLineArguments] = _

  behavior of "CromwellCommandLine"

  before {
    parser = CromwellApp.buildParser()
  }

  it should "fail to parse with no arguments" in {
    parser.parse(Array.empty[String], CommandLineArguments()).get.command shouldBe None
  }

  it should "run server when specified" in {
    parser.parse(Array("server"), CommandLineArguments()).get.command shouldBe Some(Server)
  }

  it should "fail to parse when supplying an argument to server" in {
    parser.parse(Array("server", "foo"), CommandLineArguments()) shouldBe None
  }

  it should "fail to parse with no arguments to run" in {
    parser.parse(Array("run"), CommandLineArguments()) shouldBe None
  }

  it should "fail to parse with too many arguments to run" in {
    parser.parse(Array("run", "forrest", "run"), CommandLineArguments()) shouldBe None
  }

  // --version exits the JVM which is not great in a test suite.  Haven't figure out a way to test this yet.
  //  it should "handle version output when the `-version` flag is passed" in {
  //    // I don't see a way to see that --version is printing just the version, but this at least confirms a `None`
  //    // output that should generate a usage and version.
  //    parser.parse(Array("--version"), CommandLineArguments()) shouldBe None
  //  }

  it should "run single when supplying wdl and inputs" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val optionsLast = parser.parse(Array("run", threeStep.wdl, "--inputs", threeStep.inputs), CommandLineArguments()).get
    optionsLast.command shouldBe Some(Run)
    optionsLast.workflowSource.get shouldBe threeStep.wdl
    optionsLast.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val optionsFirst = parser.parse(Array("run", "--inputs", threeStep.inputs, threeStep.wdl), CommandLineArguments()).get
    optionsFirst.command shouldBe Some(Run)
    optionsFirst.workflowSource.get shouldBe threeStep.wdl
    optionsFirst.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val validation = Try(CromwellEntryPoint.validateRunArguments(optionsFirst))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe Some(threeStep.sampleWdl.workflowSource())
    validation.get.workflowUrl shouldBe None
  }

  it should "run single when supplying workflow using url" in {
    val url = "https://path_to_url"
    val command = parser.parse(Array("run", url), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe url

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Option(url)
  }

  it should "run single when supplying workflow using url with inputs" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val url = "https://path_to_url"
    val command = parser.parse(Array("run", url, "--inputs", threeStep.inputs), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe url
    command.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Option(url)
  }

  it should "run single when supplying cwl workflow" in {
    val source = raw"""{"cwlVersion":"v1.0","class":"CommandLineTool","requirements":[{"class":"InlineJavascriptRequirement"}],"hints":[{"dockerPull":"debian:stretch-slim","class":"DockerRequirement"}],"inputs":[],"baseCommand":["touch","z","y","x","w","c","b","a"],"outputs":[{"type":"string","outputBinding":{"glob":"?","outputEval":"$${ return self.sort(function(a,b) { return a.location > b.location ? 1 : (a.location < b.location ? -1 : 0) }).map(f => f.basename).join(\" \") }\n""""
    val command = parser.parse(Array("run", "server/src/test/resources/cwl_glob_sort.cwl", "--type", "CWL"), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe "server/src/test/resources/cwl_glob_sort.cwl"

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource.get should include(source)
    validation.get.workflowUrl shouldBe None
  }

  it should "run single when supplying cwl workflow using url" in {
    val url = "https://path_to_url"
    val command = parser.parse(Array("run", url, "--type", "CWL"), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe url

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Option(url)
  }

  it should "run single when supplying wdl and inputs and options" in {
    val optionsLast = parser.parse(Array("run", "3step.wdl", "--inputs", "3step.inputs", "--options", "3step.options"), CommandLineArguments()).get
    optionsLast.command shouldBe Some(Run)
    optionsLast.workflowSource.get shouldBe "3step.wdl"
    optionsLast.workflowInputs.get.pathAsString shouldBe "3step.inputs"
    optionsLast.workflowOptions.get.pathAsString shouldBe "3step.options"

    val optionsFirst = parser.parse(Array("run", "--inputs", "3step.inputs", "--options", "3step.options", "3step.wdl"), CommandLineArguments()).get
    optionsFirst.command shouldBe Some(Run)
    optionsFirst.workflowSource.get shouldBe "3step.wdl"
    optionsFirst.workflowInputs.get.pathAsString shouldBe "3step.inputs"
    optionsFirst.workflowOptions.get.pathAsString shouldBe "3step.options"
  }

  it should "fail if workflow url is invalid" in {
    val command = parser.parse(Array("run", "htpps://url_with_invalid_protocol"), CommandLineArguments()).get
    val validation = Try(CromwellEntryPoint.validateRunArguments(command))

    validation.isFailure shouldBe true
    validation.failed.get.getMessage should include("Error while validating workflow url")
  }

  it should "fail if workflow url length is more than 2000 characters" in {
    val veryLongUrl = "https://this_url_has_more_than_2000_characters/why_would_someone_have_such_long_urls_one_would_ask/beats_me/now_starts_lorem_ipsum/At_vero_eos_et_accusamus_et_iusto_odio_dignissimos_ducimus_qui_blanditiis_praesentium_voluptatum_deleniti_atque_corrupti_quos_dolores_et_quas_molestias_excepturi_sint_occaecati_cupiditate_non_provident,_similique_sunt_in_culpa_qui_officia_deserunt_mollitia_animi,_id_est_laborum_et_dolorum_fuga._Et_harum_quidem_rerum_facilis_est_et_expedita_distinctio._Nam_libero_tempore,_cum_soluta_nobis_est_eligendi_optio_cumque_nihil_impedit_quo_minus_id_quod_maxime_placeat_facere_possimus,_omnis_voluptas_assumenda_est,_omnis_dolor_repellendus._Temporibus_autem_quibusdam_et_aut_officiis_debitis_aut_rerum_necessitatibus_saepe_eveniet_ut_et_voluptates_repudiandae_sint_et_molestiae_non_recusandae._Itaque_earum_rerum_hic_tenetur_a_sapiente_delectus,_ut_aut_reiciendis_voluptatibus_maiores_alias_consequatur_aut_perferendis_doloribus_asperiores_repellat/Sed_ut_perspiciatis_unde_omnis_iste_natus_error_sit_voluptatem_accusantium_doloremque_laudantium,_totam_rem_aperiam,_eaque_ipsa_quae_ab_illo_inventore_veritatis_et_quasi_architecto_beatae_vitae_dicta_sunt_explicabo._Nemo_enim_ipsam_voluptatem_quia_voluptas_sit_aspernatur_aut_odit_aut_fugit,_sed_quia_consequuntur_magni_dolores_eos_qui_ratione_voluptatem_sequi_nesciunt._Neque_porro_quisquam_est,_qui_dolorem_ipsum_quia_dolor_sit_amet,_consectetur,_adipisci_velit,_sed_quia_non_numquam_eius_modi_tempora_incidunt_ut_labore_et_dolore_magnam_aliquam_quaerat_voluptatem._Ut_enim_ad_minima_veniam,_quis_nostrum_exercitationem_ullam_corporis_suscipit_laboriosam,_nisi_ut_aliquid_ex_ea_commodi_consequatur?_Quis_autem_vel_eum_iure_reprehenderit_qui_in_ea_voluptate_velit_esse_quam_nihil_molestiae_consequatur,_vel_illum_qui_dolorem_eum_fugiat_quo_voluptas_nulla_pariatur?/Lorem_ipsum_dolor_sit_amet,_consectetur_adipiscing_elit,_sed_do_eiusmod_tempor_incididunt_ut_labore_et_dolore_magna_aliqua._Ut_enim_ad_minim_veniam,_quis_nostrud_exercitation_ullamco_laboris_nisi_ut_aliquip_ex_ea_commodo_consequat._Duis_aute_irure_dolor_in_reprehenderit_in_voluptate_velit_esse_cillum_dolore_eu_fugiat_nulla_pariatur._Excepteur_sint_occaecat_cupidatat_non_proident,_sunt_in_culpa_qui_officia_deserunt_mollit_anim_id_est_laborum/hello.wdl/hello.wdl"
    val command = parser.parse(Array("run", veryLongUrl), CommandLineArguments()).get
    val validation = Try(CromwellEntryPoint.validateRunArguments(command))

    validation.isFailure shouldBe true
    validation.failed.get.getMessage should include("Invalid workflow url: url has length 2305, longer than the maximum allowed 2000 characters")
  }

  it should "fail if input files do not exist" in {
    val parsedArgs = parser.parse(Array("run", "xyzshouldnotexist.wdl", "--inputs", "xyzshouldnotexist.inputs", "--options", "xyzshouldnotexist.options"), CommandLineArguments()).get
    val validation = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))

    validation.isFailure shouldBe true
    validation.failed.get.getMessage should include("Workflow source path does not exist")
    validation.failed.get.getMessage should include("Workflow inputs does not exist")
    validation.failed.get.getMessage should include("Workflow options does not exist")
  }

  it should "fail if inputs path is not readable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val parsedArgs = parser.parse(Array("run", threeStep.wdl, "--inputs", threeStep.inputs), CommandLineArguments()).get
    threeStep.inputsFile setPermissions Set.empty
    val ccl = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Workflow inputs is not readable")
  }

  it should "fail if metadata output path is not writeable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val parsedArgs = parser.parse(Array("run", threeStep.wdl, "--inputs", threeStep.inputs, "--metadata-output", threeStep.metadata), CommandLineArguments()).get
    threeStep.metadataFile write "foo"
    threeStep.metadataFile setPermissions Set.empty
    val ccl = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Unable to write to metadata directory:")
  }

  it should "run if the imports path is a .zip file" in {
    val wdlDir = DefaultPathBuilder.createTempDirectory("wdlDirectory")

    val filePassing = DefaultPathBuilder.createTempFile("filePassing", ".wdl", Option(wdlDir))
    val fileClobber = DefaultPathBuilder.createTempFile("fileClobber", ".wdl", Option(wdlDir))
    filePassing write FilePassingWorkflow.workflowSource()
    fileClobber write FileClobber.workflowSource()

    val zippedDir = wdlDir.zip()
    val zippedPath = zippedDir.pathAsString

    val parsedArgs = parser.parse(Array("run", filePassing.pathAsString, "--imports", zippedPath), CommandLineArguments()).get
    val ccl = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))
    ccl.isFailure shouldBe false

    zippedDir.delete(swallowIOExceptions = true)
  }

  it should "send content of WDL source file in Submit mode" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val command = parser.parse(Array("submit", threeStep.wdl, "--inputs", threeStep.inputs), CommandLineArguments()).get
    command.command shouldBe Some(Submit)
    command.workflowSource.get shouldBe threeStep.wdl
    command.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val validation = Try(CromwellEntryPoint.validateSubmitArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe Option(threeStep.wdlFile.contentAsString)
    validation.get.workflowUrl shouldBe None
  }
}

object CromwellCommandLineSpec {
  val ThreeStepWithoutOptions = WdlAndInputs(ThreeStep)
  val ThreeStepInputs = ThreeStepWithoutOptions.inputsFile.contentAsString

  /**
   * Create a temporary wdl file and inputs for the sampleWdl.
   * When the various properties are lazily accessed, they are also registered for deletion after the suite completes.
   */
  case class WdlAndInputs(sampleWdl: SampleWdl, optionsJson: String = "{}") {
    // Track all the temporary files we create, and delete them after the test.
    private var tempFiles = Vector.empty[Path]

    lazy val wdlFile = {
      val file = DefaultPathBuilder.createTempFile(s"${sampleWdl.name}.", ".wdl")
      tempFiles :+= file
      file write sampleWdl.workflowSource()
    }

    lazy val wdl = wdlFile.pathAsString

    lazy val inputsFile = {
      val file = wdlFile.swapExt("wdl", "inputs")
      tempFiles :+= file
      file write sampleWdl.workflowJson
    }

    lazy val inputs = inputsFile.pathAsString

    lazy val optionsFile = {
      val file = wdlFile.swapExt("wdl", "options")
      tempFiles :+= file
      file write optionsJson
    }

    lazy val options = optionsFile.pathAsString

    lazy val metadataFile = {
      val path = wdlFile.swapExt("wdl", "metadata.json")
      tempFiles :+= path
      path
    }

    lazy val metadata = metadataFile.pathAsString

    def deleteTempFiles() = tempFiles.foreach(_.delete(swallowIOExceptions = true))
  }
}
