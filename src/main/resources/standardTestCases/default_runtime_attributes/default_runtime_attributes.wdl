task default_runtime_attributes {
    command <<<
        echo "echo 'OH NO!'" > script.sh
        echo "exit 1" >> script.sh
        chmod +x script.sh
        ./script.sh
    >>>
    output {
        String ohno = read_string(stdout())
    }
    runtime {
        # Deliberately leave out the 'docker' attribute required for JES.  If that doesn't come in on the
        # default_runtime_attributes this test will fail.

        # 'continueOnReturnCode' also comes in from the default_runtime_attributes in the options file.
    }
}

workflow wf_default_runtime_attributes {
    call default_runtime_attributes
}
