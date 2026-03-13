version 1.1

# Adapted from https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#reserved-runtime-hints

task localization_tester {
  input {
    File foo
    Array[File] bar
  }

  # Require that the file is NOT localized, adapted from `draft3_nio_file_papi2.wdl`
  command <<<
    touch errors.txt
    ls

    # Check filenames are correct
    [ ~{foo} = "gs://centaur-ci/textfile_1.txt" ] || echo "Wrong filename ~{foo}" >> errors.txt
    [ ~{bar[0]} = "gs://centaur-ci/textfile_2.txt" ] || echo "Wrong filename ~{bar[0]}" >> errors.txt
    [ ~{bar[1]} = "gs://centaur-ci/textfile_3.txt" ] || echo "Wrong filename ~{bar[1]}" >> errors.txt

    # Check files are not localized
    find . -path "*centaur-ci/textfile_1.txt" >> errors.txt
    find . -path "*centaur-ci/textfile_2.txt" >> errors.txt
    find . -path "*centaur-ci/textfile_3.txt" >> errors.txt
  >>>

  output {
    String errors = read_string("errors.txt")
  }

  runtime {
    localizationOptional: true
    container: "rockylinux/rockylinux:10-minimal"
    predefinedMachineType: "e2-medium"
  }

  meta {
    volatile: true
  }
}

workflow input_hint_wf {

  input {
    File foo
    Array[File] bar
  }

  call localization_tester {
    input: foo, bar
  }

  output {
    String errors = localization_tester.errors
  }
}
