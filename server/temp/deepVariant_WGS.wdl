version 1.0

task MakeExamples {
    input {
        File input_bam
        File input_bam_index

        Int startIndex
        Int endIndex
        Int shards
        String? regions

        File input_ref
        File input_ref_fai
        File? input_ref_gzi

        String docker_image
        Float ram
        Int cores
        Int preemptible_tries
    }

    Int disk_size = ceil((size(input_bam, "GB") + size(input_ref, "GB"))*2 + 50)

    command <<<
        mkdir examples buffers

        seq "~{startIndex}" "~{endIndex}" | parallel --halt 2 \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --examples "examples/examples_output.tfrecord@~{shards}".gz \
        --gvcf "buffers/gvcf_output.tfrecord@~{shards}".gz \
        --reads "~{input_bam}" \
        --ref "~{input_ref}" \
        --task {} ~{"--regions " + regions}

        tar -zcvf buffers.tar.gz buffers
    >>>

    output {
        Array[File] examples = glob("examples/*")
        File variantProtocolBuffers = "buffers.tar.gz"
    }

    runtime {
        docker: docker_image
        cpu: cores
        memory: ram + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
    }
}

task WriteListOfFiles {
    input {
        Array[Array[String]] input_files
        File tsv = write_tsv(input_files)

        Int preemptible_tries
    }

    command <<<
        python <<CODE
        with open("~{tsv}","r") as tsv_file:
            for line in tsv_file:
                parts = line.strip().split("\t")
                for part in parts:
                    if (part[-2:] == "gz"):
                        print part
        CODE
    >>>

    output {
        File list_out = stdout()
    }

    runtime {
        docker: "python:2.7"
        cpu: 1
        memory: "0.6 GB"
        preemptible: preemptible_tries
    }
}

task GetExamples {
    input {
        File examples
        Int start_index
        Int end_index

        Int preemptible_tries
    }

    command <<<
        python <<CODE
        with open("~{examples}","r") as examples:
            start_index = ~{start_index}
            end_index = ~{end_index}
            i = 0
            for line in examples:
                if i >= start_index and i <= end_index:
                    print line.strip()
                i += 1

        CODE
    >>>

    output {
        Array[String] examples_for_shard = read_lines(stdout())
    }

    runtime {
        docker: "python:2.7"
        cpu: 1
        memory: "0.6 GB"
        preemptible: preemptible_tries
    }
}

task CallVariants {
    input {
        Array[File] examples
        Int startIndex
        Int endIndex
        Int concurrent_jobs
        Int shards

        File model
        File model_index
        File model_meta

        String docker_image
        Int cores
        Int ram
        Int preemptible_tries
    }

    command <<<
        mkdir called_variants
        mkdir examples
        mkdir model

        for i in $(echo ~{sep=" "examples}); do
        cp  $i examples/
        done

        MODEL_DIR=$(dirname ~{model})

        seq -f "%05g" "~{startIndex}" "~{endIndex}" | \
        parallel --jobs "~{concurrent_jobs}" --halt 2 \
        /opt/deepvariant/bin/call_variants \
        --examples examples/examples_output.tfrecord-{}-of-"$(printf "%05d" "~{shards}")".gz \
        --outfile called_variants/call_variants_output.tfrecord-{}-of-"$(printf "%05d" "~{shards}")".gz \
        --checkpoint $MODEL_DIR/model.ckpt

        tar -zcvf called_variants.tar.gz called_variants
    >>>

    output {
        File called_variants_tar = "called_variants.tar.gz"
    }

    runtime {
        docker: docker_image
        cpu: cores
        memory: ram + " GB"
        disks: "local-disk 75 HDD"
        preemptible: preemptible_tries
    }
}

task gatherFiles {
    input {
        Array [File] tarfiles
        String docker_image
        Int preemptible_tries
    }

    String base = basename(tarfiles[0], ".tar.gz")

    command {
        for tarfile in $(echo -e "~{sep='\n' tarfiles}"); do
        tar -zxvf $tarfile
        done

        tar -zcvf ~{base}.tar.gz ~{base}
    }

    output {
        File gathered = "~{base}.tar.gz"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk 50 HDD"
        preemptible: preemptible_tries
    }
}

task PostProcessVariants {
    input {
        File called_variants_all
        File variantProtocolBuffers
        String out_vcf_file_name
        String out_gvcf_file_name

        File input_ref
        File input_ref_index
        File? input_ref_gzi

        String docker_image
        Int shards
        Int preemptible_tries
    }

    command <<<
        tar -zxvf ~{called_variants_all}
        tar -zxvf ~{variantProtocolBuffers}

        /opt/deepvariant/bin/postprocess_variants \
        --ref "~{input_ref}" \
        --infile "called_variants/call_variants_output.tfrecord@~{shards}".gz \
        --nonvariant_site_tfrecord_path "buffers/gvcf_output.tfrecord@~{shards}".gz \
        --gvcf_outfile "~{out_gvcf_file_name}" \
        --outfile "~{out_vcf_file_name}"
    >>>

    output {
        File output_vcf = "~{out_vcf_file_name}"
        File vcf_index = "~{out_vcf_file_name}.tbi"
        File output_gvcf = "~{out_gvcf_file_name}"
        File gvcf_index = "~{out_gvcf_file_name}.tbi"
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "30 GB"
        disks: "local-disk 50 HDD"
        preemptible: preemptible_tries
    }
}

workflow WGSDeepVariantCostOpt {
    input {
        String input_bam
        String input_bam_index

        String input_ref
        String input_ref_index
        String? input_ref_gzi

        String? regions
    }

    File model = "models_DeepVariant_0.10.0_DeepVariant-inception_v3-0.10.0+data-wgs_standard_model.ckpt.data-00000-of-00001"
    File model_index = "model.ckpt.index"
    File model_meta = "model.ckpt.meta"

    String image_version = "0.10.0"
    String docker_image = "gcr.io/deepvariant-docker/deepvariant:" + image_version
    String dnastack_project = "cool-benefit-817"
    Int preemptible_tries = 5

    String output_basename = basename(input_bam, ".bam")

    Int shards = 128
    Int num_cores_per_shard = 4

    Int num_example_workers = 8
    Int num_shards_per_worker_example = shards / num_example_workers

    Int example_cores = 16
    Float example_ram = 36.0

    Int num_variants_workers = 8
    Int num_shards_per_worker_variants = shards / num_variants_workers

    Int variants_cores = 32
    Int variants_ram = 60

    Int concurrent_jobs = variants_cores / num_cores_per_shard

    scatter(worker in range(num_example_workers)) {
        Int example_start_index = worker * num_shards_per_worker_example
        Int example_end_index = ((worker + 1) * num_shards_per_worker_example ) - 1

        call MakeExamples{
            input:
                startIndex = example_start_index,
                endIndex = example_end_index,
                input_bam = input_bam,
                input_bam_index = input_bam_index,
                input_ref = input_ref,
                input_ref_fai = input_ref_index,
                shards = shards,
                cores = example_cores,
                ram = example_ram,
                docker_image = docker_image,
                preemptible_tries = preemptible_tries,
                regions = regions,
                input_ref_gzi = input_ref_gzi
        }
    }

    call WriteListOfFiles {
        input:
            input_files = MakeExamples.examples,
            preemptible_tries = preemptible_tries
    }

    scatter (worker in range(num_variants_workers)) {
        Int variants_start_index = worker * num_shards_per_worker_variants
        Int variants_end_index = ((worker + 1) * num_shards_per_worker_variants) - 1

        call GetExamples {
            input:
                examples = WriteListOfFiles.list_out,
                start_index = variants_start_index,
                end_index = variants_end_index,
                preemptible_tries = preemptible_tries
        }

        call CallVariants {
            input:
                startIndex = variants_start_index,
                endIndex = variants_end_index,
                examples = GetExamples.examples_for_shard,
                model = model,
                model_index = model_index,
                model_meta = model_meta,
                shards = shards,
                concurrent_jobs = concurrent_jobs,
                cores = variants_cores,
                ram = variants_ram,
                docker_image = docker_image,
                preemptible_tries = preemptible_tries
        }
    }

    call gatherFiles as gatherVariants {
        input:
            tarfiles = CallVariants.called_variants_tar,
            preemptible_tries = preemptible_tries,
            docker_image = docker_image
    }

    call gatherFiles as gatherBuffers {
        input:
            tarfiles = MakeExamples.variantProtocolBuffers,
            preemptible_tries = preemptible_tries,
            docker_image = docker_image
    }

    call PostProcessVariants {
        input:
            called_variants_all = gatherVariants.gathered,
            variantProtocolBuffers = gatherBuffers.gathered,
            input_ref = input_ref,
            input_ref_index = input_ref_index,
            shards = shards,
            out_vcf_file_name = output_basename + ".deep_variant.vcf.gz",
            out_gvcf_file_name = output_basename + ".deep_variant.g.vcf.gz",
            docker_image = docker_image,
            preemptible_tries = preemptible_tries,
            input_ref_gzi = input_ref_gzi
    }


    output {
        File output_vcf = PostProcessVariants.output_vcf
        File output_gvcf = PostProcessVariants.output_gvcf
        File vcf_index = PostProcessVariants.vcf_index
        File gvcf_index = PostProcessVariants.gvcf_index
    }

    meta {
        author: "Patrick Magee"
        email: "patrick@dnastack.com"
        version: "0005"
        last_modified: "2020-03-30"
        last_revised_by: "Heather Ward"
    }
}