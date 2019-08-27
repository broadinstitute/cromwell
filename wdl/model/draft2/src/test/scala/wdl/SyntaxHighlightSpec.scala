package wdl

import org.scalatest.{Matchers, WordSpec}
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model.WdlNamespace
import wdl.draft2.model.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter}
import wom.ResolvedImportRecord

class SyntaxHighlightSpec extends WordSpec with Matchers {
  "SyntaxFormatter for typical workflow" should {
    val namespace = WdlNamespace.loadUsingSource(
      """
        |task PairedFastQsToUnmappedBAM {
        |  File fastq_1
        |  File fastq_2
        |  String readgroup_name
        |  String sample_name
        |  String library_name
        |  String platform_unit
        |  String run_date
        |  String platform_name
        |  String sequencing_center
        |  Int disk_size
        |  String mem_size
        |
        |  command {
        |    java -Xmx3000m -jar /usr/gitc/picard.jar \
        |      FastqToSam \
        |      FASTQ=${fastq_1} \
        |      FASTQ2=${fastq_2} \
        |      OUTPUT=${readgroup_name}.bam \
        |      READ_GROUP_NAME=${readgroup_name} \
        |      SAMPLE_NAME=${sample_name} \
        |      LIBRARY_NAME=${library_name} \
        |      PLATFORM_UNIT=${platform_unit} \
        |      RUN_DATE=${run_date} \
        |      PLATFORM=${platform_name} \
        |      SEQUENCING_CENTER=${sequencing_center}
        |  }
        |  runtime {
        |    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        |    memory: mem_size
        |    cpu: "1"
        |    disks: "local-disk " + disk_size + " HDD"
        |  }
        |  output {
        |    File output_bam = "${readgroup_name}.bam"
        |  }
        |  	parameter_meta {
        |    memory_mb: "Amount of memory to allocate to the JVM"
        |    param: "Some arbitrary parameter"
        |    sample_id: "The ID of the sample in format foo_bar_baz"
        |  }
        |  meta {
        |    author: "Joe Somebody"
        |    email: "joe@company.org"
        |  }
        |}
        |
        |# WORKFLOW DEFINITION
        |workflow ConvertPairedFastQsToUnmappedBamWf {
        |  Array[String] readgroup_list
        |  Map[String, Array[File]] fastq_pairs
        |  Map[String, Array[String]] metadata
        |
        |  # Convert multiple pairs of input fastqs in parallel
        |  scatter (readgroup in readgroup_list) {
        |
        |    # Convert pair of FASTQs to uBAM
        |    call PairedFastQsToUnmappedBAM {
        |      input:
        |        fastq_1 = fastq_pairs[readgroup][0],
        |        fastq_2 = fastq_pairs[readgroup][1],
        |        readgroup_name = readgroup,
        |        sample_name = metadata[readgroup][0],
        |        library_name = metadata[readgroup][1],
        |        platform_unit = metadata[readgroup][2],
        |        run_date = metadata[readgroup][3],
        |        platform_name = metadata[readgroup][4],
        |        sequencing_center = metadata[readgroup][5]
        |    }
        |  }
        |
        |  # Outputs that will be retained when execution is complete
        |  output {
        |    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
        |  }
        |
        |  	parameter_meta {
        |    memory_mb: "Amount of memory to allocate to the JVM"
        |    param: "Some arbitrary parameter"
        |    sample_id: "The ID of the sample in format foo_bar_baz"
        |  }
        |
        |  meta {
        |   author: "Joe Somebody"
        |   email: "joe@company.org"
        |  }
        |}
     """.stripMargin, None, None).get

    val console =
      """\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mPairedFastQsToUnmappedBAM\u001b[0m {
        |  \u001b[38;5;33mFile\u001b[0m \u001b[38;5;112mfastq_1\u001b[0m
        |  \u001b[38;5;33mFile\u001b[0m \u001b[38;5;112mfastq_2\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mreadgroup_name\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112msample_name\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mlibrary_name\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mplatform_unit\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mrun_date\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mplatform_name\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112msequencing_center\u001b[0m
        |  \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mdisk_size\u001b[0m
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mmem_size\u001b[0m
        |  \u001b[38;5;214mcommand\u001b[0m {
        |    java -Xmx3000m -jar /usr/gitc/picard.jar \
        |      FastqToSam \
        |      FASTQ=${fastq_1} \
        |      FASTQ2=${fastq_2} \
        |      OUTPUT=${readgroup_name}.bam \
        |      READ_GROUP_NAME=${readgroup_name} \
        |      SAMPLE_NAME=${sample_name} \
        |      LIBRARY_NAME=${library_name} \
        |      PLATFORM_UNIT=${platform_unit} \
        |      RUN_DATE=${run_date} \
        |      PLATFORM=${platform_name} \
        |      SEQUENCING_CENTER=${sequencing_center}
        |  }
        |  \u001b[38;5;214moutput\u001b[0m {
        |    \u001b[38;5;33mFile\u001b[0m \u001b[38;5;112moutput_bam\u001b[0m = "${readgroup_name}.bam"
        |  }
        |  \u001b[38;5;214mruntime\u001b[0m {
        |    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        |    memory: mem_size
        |    cpu: "1"
        |    disks: "local-disk " + disk_size + " HDD"
        |  }
        |  \u001b[38;5;214mmeta\u001b[0m {
        |    author: "Joe Somebody"
        |    email: "joe@company.org"
        |  }
        |  \u001b[38;5;214mparameter_meta\u001b[0m {
        |    memory_mb: "Amount of memory to allocate to the JVM"
        |    param: "Some arbitrary parameter"
        |    sample_id: "The ID of the sample in format foo_bar_baz"
        |  }
        |}
        |
        |\u001b[38;5;214mworkflow\u001b[0m \u001b[38;5;253mConvertPairedFastQsToUnmappedBamWf\u001b[0m {
        |  \u001b[38;5;33mArray[String]\u001b[0m \u001b[38;5;112mreadgroup_list\u001b[0m
        |  \u001b[38;5;33mMap[String, Array[File]]\u001b[0m \u001b[38;5;112mfastq_pairs\u001b[0m
        |  \u001b[38;5;33mMap[String, Array[String]]\u001b[0m \u001b[38;5;112mmetadata\u001b[0m
        |  \u001b[38;5;214mscatter\u001b[0m (readgroup in readgroup_list) {
        |    \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mPairedFastQsToUnmappedBAM\u001b[0m {
        |      input: library_name=metadata[readgroup][1], run_date=metadata[readgroup][3], readgroup_name=readgroup, platform_name=metadata[readgroup][4], platform_unit=metadata[readgroup][2], fastq_1=fastq_pairs[readgroup][0], fastq_2=fastq_pairs[readgroup][1], sample_name=metadata[readgroup][0], sequencing_center=metadata[readgroup][5]
        |    }
        |  }
        |  \u001b[38;5;33mArray[File]\u001b[0m \u001b[38;5;112moutput_bams\u001b[0m = PairedFastQsToUnmappedBAM.output_bam
        |  \u001b[38;5;214mmeta\u001b[0m {
        |    author: "Joe Somebody"
        |    email: "joe@company.org"
        |  }
        |  \u001b[38;5;214mparameter_meta\u001b[0m {
        |    memory_mb: "Amount of memory to allocate to the JVM"
        |    param: "Some arbitrary parameter"
        |    sample_id: "The ID of the sample in format foo_bar_baz"
        |  }
        |}""".stripMargin

    val html =
      """<span class="keyword">task</span> <span class="name">PairedFastQsToUnmappedBAM</span> {
        |  <span class="type">File</span> <span class="variable">fastq_1</span>
        |  <span class="type">File</span> <span class="variable">fastq_2</span>
        |  <span class="type">String</span> <span class="variable">readgroup_name</span>
        |  <span class="type">String</span> <span class="variable">sample_name</span>
        |  <span class="type">String</span> <span class="variable">library_name</span>
        |  <span class="type">String</span> <span class="variable">platform_unit</span>
        |  <span class="type">String</span> <span class="variable">run_date</span>
        |  <span class="type">String</span> <span class="variable">platform_name</span>
        |  <span class="type">String</span> <span class="variable">sequencing_center</span>
        |  <span class="type">Int</span> <span class="variable">disk_size</span>
        |  <span class="type">String</span> <span class="variable">mem_size</span>
        |  <span class="section">command</span> {
        |    <span class="command">java -Xmx3000m -jar /usr/gitc/picard.jar \
        |      FastqToSam \
        |      FASTQ=${fastq_1} \
        |      FASTQ2=${fastq_2} \
        |      OUTPUT=${readgroup_name}.bam \
        |      READ_GROUP_NAME=${readgroup_name} \
        |      SAMPLE_NAME=${sample_name} \
        |      LIBRARY_NAME=${library_name} \
        |      PLATFORM_UNIT=${platform_unit} \
        |      RUN_DATE=${run_date} \
        |      PLATFORM=${platform_name} \
        |      SEQUENCING_CENTER=${sequencing_center}</span>
        |  }
        |  <span class="section">output</span> {
        |    <span class="type">File</span> <span class="variable">output_bam</span> = "${readgroup_name}.bam"
        |  }
        |  <span class="keyword">runtime</span> {
        |    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        |    memory: mem_size
        |    cpu: "1"
        |    disks: "local-disk " + disk_size + " HDD"
        |  }
        |  <span class="keyword">meta</span> {
        |    author: "Joe Somebody"
        |    email: "joe@company.org"
        |  }
        |  <span class="keyword">parameter_meta</span> {
        |    memory_mb: "Amount of memory to allocate to the JVM"
        |    param: "Some arbitrary parameter"
        |    sample_id: "The ID of the sample in format foo_bar_baz"
        |  }
        |}
        |
        |<span class="keyword">workflow</span> <span class="name">ConvertPairedFastQsToUnmappedBamWf</span> {
        |  <span class="type">Array[String]</span> <span class="variable">readgroup_list</span>
        |  <span class="type">Map[String, Array[File]]</span> <span class="variable">fastq_pairs</span>
        |  <span class="type">Map[String, Array[String]]</span> <span class="variable">metadata</span>
        |  <span class="keyword">scatter</span> (readgroup in readgroup_list) {
        |    <span class="keyword">call</span> <span class="name">PairedFastQsToUnmappedBAM</span> {
        |      input: library_name=metadata[readgroup][1], run_date=metadata[readgroup][3], readgroup_name=readgroup, platform_name=metadata[readgroup][4], platform_unit=metadata[readgroup][2], fastq_1=fastq_pairs[readgroup][0], fastq_2=fastq_pairs[readgroup][1], sample_name=metadata[readgroup][0], sequencing_center=metadata[readgroup][5]
        |    }
        |  }
        |  <span class="type">Array[File]</span> <span class="variable">output_bams</span> = PairedFastQsToUnmappedBAM.output_bam
        |  <span class="keyword">meta</span> {
        |    author: "Joe Somebody"
        |    email: "joe@company.org"
        |  }
        |  <span class="keyword">parameter_meta</span> {
        |    memory_mb: "Amount of memory to allocate to the JVM"
        |    param: "Some arbitrary parameter"
        |    sample_id: "The ID of the sample in format foo_bar_baz"
        |  }
        |}""".stripMargin

    "format to console properly" in {
      val actual = new SyntaxFormatter(AnsiSyntaxHighlighter).format(namespace)
      actual shouldEqual console
    }

    "format to HTML properly" in {
      val actual = new SyntaxFormatter(HtmlSyntaxHighlighter).format(namespace)
      actual shouldEqual html
    }
  }

  "SyntaxFormatter for more feature-rich workflow" should {

    val fooTaskWdl = """
                       |task foo {
                       |  command {
                       |    echo "foo!"
                       |  }
                       |  output {
                       |  File out = stdout()
                       |  }
                       |}
                     """.stripMargin

    def resolver(importUri: String): Draft2ResolvedImportBundle = {
      importUri match {
        case "foo.wdl" => Draft2ResolvedImportBundle(fooTaskWdl, ResolvedImportRecord("foo.wdl"))
        case _ => throw new RuntimeException(s"Can't resolve $importUri")
      }
    }

    val source = s"""
        |import "foo.wdl" as foo_ns
        |
        |task t {
        |  String f
        |  Int p
        |  command {
        |    ./cmd $${f} $${p}
        |  }
        |}
        |
        |task s {
        |  Array[File] input_file
        |  command <<<
        |    cat $${sep=' ' input_file} | awk '{s+=$$1} END {print s}'
        |  >>>
        |  output {
        |    String s = read_string(stdout())
        |  }
        |}
        |
        |task r {
        |  command { python -c "import random; print(random.randint(1,100))" }
        |}
        |
        |workflow w {
        |  Int p = 2+2
        |  call t
        |  call t as u {
        |    input: f="abc", p=p
        |  }
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(source, None, Option(Seq(resolver))).get

    val console =
      s"""\u001b[38;5;214mimport\u001b[0m 'foo.wdl' as foo_ns
        |
        |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mt\u001b[0m {
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mf\u001b[0m
        |  \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mp\u001b[0m
        |  \u001b[38;5;214mcommand\u001b[0m {
        |    ./cmd $${f} $${p}
        |  }
        |}
        |
        |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253ms\u001b[0m {
        |  \u001b[38;5;33mArray[File]\u001b[0m \u001b[38;5;112minput_file\u001b[0m
        |  \u001b[38;5;214mcommand\u001b[0m <<<
        |    cat $${sep=" " input_file} | awk '{s+=$$1} END {print s}'
        |  >>>
        |  \u001b[38;5;214moutput\u001b[0m {
        |    \u001b[38;5;33mString\u001b[0m \u001b[38;5;112ms\u001b[0m = \u001b[38;5;13mread_string\u001b[0m(\u001b[38;5;13mstdout\u001b[0m())
        |  }
        |}
        |
        |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mr\u001b[0m {
        |  \u001b[38;5;214mcommand\u001b[0m {
        |    python -c "import random; print(random.randint(1,100))"
        |  }
        |}
        |
        |\u001b[38;5;214mworkflow\u001b[0m \u001b[38;5;253mw\u001b[0m {
        |  \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mp\u001b[0m = 2 + 2
        |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m
        |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m as u {
        |    input: f="abc", p=p
        |  }
        |}""".stripMargin

    val html =
      s"""<span class="keyword">import</span> 'foo.wdl' as foo_ns
        |
        |<span class="keyword">task</span> <span class="name">t</span> {
        |  <span class="type">String</span> <span class="variable">f</span>
        |  <span class="type">Int</span> <span class="variable">p</span>
        |  <span class="section">command</span> {
        |    <span class="command">./cmd $${f} $${p}</span>
        |  }
        |}
        |
        |<span class="keyword">task</span> <span class="name">s</span> {
        |  <span class="type">Array[File]</span> <span class="variable">input_file</span>
        |  <span class="section">command</span> <<<
        |    <span class="command">cat $${sep=" " input_file} | awk '{s+=$$1} END {print s}'</span>
        |  >>>
        |  <span class="section">output</span> {
        |    <span class="type">String</span> <span class="variable">s</span> = <span class="function">read_string</span>(<span class="function">stdout</span>())
        |  }
        |}
        |
        |<span class="keyword">task</span> <span class="name">r</span> {
        |  <span class="section">command</span> {
        |    <span class="command">python -c "import random; print(random.randint(1,100))"</span>
        |  }
        |}
        |
        |<span class="keyword">workflow</span> <span class="name">w</span> {
        |  <span class="type">Int</span> <span class="variable">p</span> = 2 + 2
        |  <span class="keyword">call</span> <span class="name">t</span>
        |  <span class="keyword">call</span> <span class="name">t</span> as <span class="alias">u</span> {
        |    input: f="abc", p=p
        |  }
        |}""".stripMargin

    "format to console properly" in {
      val actual = new SyntaxFormatter(AnsiSyntaxHighlighter).format(namespace)
      actual shouldEqual console
    }

    "format to HTML properly" in {
      val actual = new SyntaxFormatter(HtmlSyntaxHighlighter).format(namespace)
      actual shouldEqual html
    }
  }
}
