version 1.0

workflow cwl_conformance_test {
    input {
        String cwl_dir
        String test_result_output
        String centaur_cwl_runner
        String conformance_expected_failures
        Int timeout
    }

    call get_test_count {
        input:
            cwl_dir = cwl_dir,
            centaur_cwl_runner = centaur_cwl_runner
    }

    scatter (test_index in range(get_test_count.test_count)) {
        call run_test_index {
            input:
                cwl_dir = cwl_dir,
                centaur_cwl_runner = centaur_cwl_runner,
                test_index = test_index,
                timeout = timeout
        }
    }

    call make_summary {
        input:
            test_count = get_test_count.test_count,
            test_result_codes = run_test_index.test_result_code,
            conformance_expected_failures = conformance_expected_failures,
            test_result_outputs = run_test_index.out,
            test_result_output = test_result_output
    }

    output {
    }
}

task get_test_count {
    input {
        String cwl_dir
        String centaur_cwl_runner
    }

    command {
        cd ~{cwl_dir}
        ./run_test.sh RUNNER="~{centaur_cwl_runner}" -l | grep -c '^\['
    }

    output {
        Int test_count = read_int(stdout())
    }
}

task run_test_index {
    input {
        String cwl_dir
        String centaur_cwl_runner
        Int test_index
        Int test_number = test_index + 1
        Int timeout
    }

    # Weird -n/--timeout format is because ./run_test.sh doesn't pass through the timeout parameter to cwltest.
    # Test 55 often runs over 10 minutes in Travis/PapiV2 and requires a longer timeout.
    # So, we use the fact that run_test.sh conveniently doesn't sanitize -n to wire the --timeout _inside_ the -n arg.
    command {
        (
            cd ~{cwl_dir}
            ./run_test.sh RUNNER="~{centaur_cwl_runner}" -n"~{test_number} --timeout=~{timeout}" 2>&1
        )
        echo $? > test_result_code
    }

    output {
        Int test_result_code = read_int("test_result_code")
        File out = stdout()
    }
}

task make_summary {
    input {
        Int test_count
        Array[String] test_result_outputs
        Array[Int] test_result_codes
        String test_result_output
        String conformance_expected_failures
        File test_result_output_lines = write_lines(test_result_outputs)
        File test_result_code_lines = write_lines(test_result_codes)
    }

    command <<<
        TEST_PASSING=0
        touch unexpected_pass
        touch unexpected_fail

        for TEST_NUMBER in $(seq ~{test_count}); do
            # Check if test is supposed to fail
            grep -q "^${TEST_NUMBER}$" "~{conformance_expected_failures}"
            TEST_IN_EXPECTED_FAILED=$?

            # Get the test results
            TEST_RESULT_OUTPUT="$(sed -n ${TEST_NUMBER}p ~{test_result_output_lines})"
            TEST_RESULT_CODE="$(sed -n ${TEST_NUMBER}p ~{test_result_code_lines})"

            if [ "${TEST_RESULT_CODE}" -eq 0 ]; then
                TEST_PASSING="$((${TEST_PASSING} + 1))"
            fi

            # Check for unexpected results
            if [ "${TEST_IN_EXPECTED_FAILED}" -eq 0 ] && [ "${TEST_RESULT_CODE}" -eq 0 ]; then
                echo "${TEST_NUMBER}" >> unexpected_pass
            elif [ "${TEST_IN_EXPECTED_FAILED}" -ne 0 ] && [ "${TEST_RESULT_CODE}" -ne 0 ]; then
                echo "${TEST_NUMBER}" >> unexpected_fail
            fi

            cat "$TEST_RESULT_OUTPUT" >> "~{test_result_output}"
            echo "exited with code ${TEST_RESULT_CODE}" >> "~{test_result_output}"
        done

        echo "---" >> "~{test_result_output}"
        echo "Conformance percentage at $(( 100 * ${TEST_PASSING} / ~{test_count} ))%" >> "~{test_result_output}"

        if [ -s unexpected_pass ] || [ -s unexpected_fail ]; then
            printf "Unexpected passing tests: (%s)\n" "$(paste -s -d ' ' unexpected_pass)" >> "~{test_result_output}"
            printf "Unexpected failing tests: (%s)\n" "$(paste -s -d ' ' unexpected_fail)" >> "~{test_result_output}"
            echo "Does ~{conformance_expected_failures} need to be updated?" >> "~{test_result_output}"
            false
        fi
    >>>

    output {
    }
}
