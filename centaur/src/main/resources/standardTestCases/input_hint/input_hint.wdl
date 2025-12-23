version development-1.1

task test_hints {
  input {
    File foo
    File bar
    File baz
    File qux
  }

  # The version of this workflow in the WDL spec does a `wc` on a localized file.
  # I guess it demonstrates an engine that ignores `localizationOptional`.
  #
  # We are testing the opposite, require that the file is NOT localized.
  # Adapted from `draft3_nio_file_papi2.wdl` and
  # `https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#reserved-runtime-hints`
  command <<<
    echo ~{foo} >> filename1.txt
    find . -name ~{foo} >> errors.txt

    echo ~{bar} >> filename2.txt
    echo ~{baz} >> filename3.txt
    echo ~{qux} >> filename4.txt
  >>>

  output {
    String filename1 = read_string("filename1.txt")
    String filename2 = read_string("filename2.txt")
    String filename3 = read_string("filename3.txt")
    String filename4 = read_string("filename4.txt")
    String errors = read_string("errors.txt")
  }

  runtime {
    container: "ubuntu:latest"
    maxMemory: "36 GB"
    maxCpu: 24
    shortTask: true
    localizationOptional: false
    inputs: object {
      foo: object {
        localizationOptional: true
      },
      bar: object {
        localizationOptional: false
      },
      baz: object {
        banana: true
      },
      qux: object {
        banana: "Zardoz"
      }
    }
  }

  meta {
    volatile: true
  }
}

workflow input_hint_wf {

  call test_hints {
    input: foo = "gs://centaur-ci/textfile_1.txt",
           bar = "gs://centaur-ci/textfile_2.txt",
           baz = "gs://centaur-ci/textfile_3.txt",
           qux = "gs://centaur-ci/textfile_4.txt"
  }

  output {
    String errors = test_hints.errors
    String filename1 = test_hints.filename1
    String filename2 = test_hints.filename2
    String filename3 = test_hints.filename3
    String filename4 = test_hints.filename4
  }
}
