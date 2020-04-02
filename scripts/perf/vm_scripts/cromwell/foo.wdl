version 1.0
workflow build_augur_tree {
	input {
		# fasta header records need to be "|" delimited for each metadata value
		File raw_sequences
		File ref_fasta
		# list metadata values space delimited found in headers of raw_sequences.fasta
		String fasta_fields
	}
	call metadata_parse {
		input:
			raw_sequences = raw_sequences,
			fasta_fields = fasta_fields
	}
    call align_sequences {
		input:
			parsed_sequences = metadata_parse.parsed_sequences_fasta,
			ref_fasta = ref_fasta
	}
    call build_tree {
		input:
			aligned_fasta = align_sequences.aligned_sequences
	}
    #output {
    #	File metadata_tsv = metadata_parse.parsed_metadata_tsv
	#	File sequences_only_fasta = metadata_parse.parsed_sequences_fasta
    #    File aligned_fasta = align_sequences.aligned_sequences
    #}
}
task metadata_parse {
	input {
		File raw_sequences
		String fasta_fields
	}
    	String raw_sequences_basename = basename(raw_sequences, ".fasta")
	command {
		augur parse --sequences ~{raw_sequences} \
			--output-sequences ~{raw_sequences_basename}_parsed_sequences.fasta \
			--output-metadata ~{raw_sequences_basename}_parsed_metadata.tsv \
			--fields ~{fasta_fields}
	}
	runtime {
		docker: "nextstrain/base:latest"
		memory: "5GB"
        preemptible: 3
	}
	output {
		File parsed_metadata_tsv = "${raw_sequences_basename}_parsed_metadata.tsv"
		File parsed_sequences_fasta = "${raw_sequences_basename}_parsed_sequences.fasta"
	}
}
task align_sequences {
	input {
		File parsed_sequences
		File ref_fasta
	}
	command {
		augur align --sequences ~{parsed_sequences} \
			--reference-sequence ~{ref_fasta} \
			--output ~{parsed_sequences}_aligned.fasta \
			--fill-gaps \
			--remove-reference \
			--nthreads auto
	}
	runtime {
		docker: "nextstrain/base:latest"
		memory: "5GB"
		preemptible: 3
	}
	output {
		File aligned_sequences = "~{parsed_sequences}_aligned.fasta"
	}
}
task build_tree {
	input {
		File aligned_fasta
		String raw_sequences_basename
	}
	command {
		augur tree --alignment ~{aligned_fasta} \
			--output ~{raw_sequences_basename}.tree \
			--nthreads auto
	}
	runtime {
		docker: "nextstrain/base:latest"
		memory: "5GB"
		preemptible: 3
	}
	output {
		File aigned_tree = "~{raw_sequences_basename}.tree"
	}
}
