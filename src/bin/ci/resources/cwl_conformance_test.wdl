workflow cwl_conformance_test {
    String cwl_dir
    String test_result_output
    String centaur_cwl_runner
    String conformance_expected_failures

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
                test_index = test_index
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

    call echo_summary {
        input:
            summary_result_text = make_summary.summary_result_text,
            summary_result_code = make_summary.summary_result_code
    }

    output {
        String summary_result_text = make_summary.summary_result_text
        Int summary_result_code = make_summary.summary_result_code
    }
}

task get_test_count {
    String cwl_dir
    String centaur_cwl_runner

    command {
        cd ${cwl_dir}
        ./run_test.sh RUNNER="${centaur_cwl_runner}" -l | grep '^\[' | wc -l
    }

    output {
        Int test_count = read_int(stdout())
    }
}

task run_test_index {
    String cwl_dir
    String centaur_cwl_runner
    Int test_index
    Int test_number = test_index + 1

    command {
        (
            cd ${cwl_dir}
            ./run_test.sh RUNNER="${centaur_cwl_runner}" -n${test_number} 2>&1
        )
        echo $? > test_result_code
    }

    output {
        Int test_result_code = read_int("test_result_code")
        File out = stdout()
    }
}

task make_summary {
    Int test_count
    Array[String] test_result_outputs
    Array[Int] test_result_codes
    String test_result_output
    String conformance_expected_failures
    File test_result_output_lines = write_lines(test_result_outputs)
    File test_result_code_lines = write_lines(test_result_codes)
    String varBegin = "${"
    String varEnd = "}"

    command {
        TEST_PASSING=0
        touch unexpected_pass
        touch unexpected_fail
        for TEST_NUMBER in $(seq ${test_count}); do
            # Check if test is supposed to fail
            grep -q '^'$TEST_NUMBER'$' ${conformance_expected_failures}
            TEST_IN_EXPECTED_FAILED=$?

            # Get the test results
            TEST_RESULT_OUTPUT=$(sed -n ${varBegin}TEST_NUMBER${varEnd}p ${test_result_output_lines})
            TEST_RESULT_CODE=$(sed -n ${varBegin}TEST_NUMBER${varEnd}p ${test_result_code_lines})

            if [ $TEST_RESULT_CODE -eq 0 ]; then
                TEST_PASSING=$(($TEST_PASSING + 1))
            fi

            # Check for unexpected results
            if [ $TEST_IN_EXPECTED_FAILED -eq 0 ] && [ $TEST_RESULT_CODE -eq 0 ]; then
                echo $TEST_NUMBER >> unexpected_pass
            elif [ ! $TEST_IN_EXPECTED_FAILED -eq 0 ] && [ ! $TEST_RESULT_CODE -eq 0 ]; then
                echo $TEST_NUMBER >> unexpected_fail
            fi

            cat "$TEST_RESULT_OUTPUT" >> "${test_result_output}"
            echo exited with code "$TEST_RESULT_CODE" >> "${test_result_output}"
        done
        echo Conformance percentage at $(( 100 * $TEST_PASSING / ${test_count} ))% >> summary_result_text

        if [ -s unexpected_pass ] || [ -s unexpected_fail ]; then
            printf 'Unexpected passing tests: (%s)\n' "$(paste -s -d ' ' unexpected_pass)" >> summary_result_text
            printf 'Unexpected failing tests: (%s)\n' "$(paste -s -d ' ' unexpected_fail)" >> summary_result_text
            printf 'Does ${conformance_expected_failures} need to be updated?\n' >> summary_result_text
            echo 1 > summary_result_code
        else
            echo 0 > summary_result_code
        fi
    }

    output {
        String summary_result_text = read_string("summary_result_text")
        Int summary_result_code = read_int("summary_result_code")
    }
}

task echo_summary {
    String summary_result_text
    Int summary_result_code

    command {
        echo <<SUMMARY
        ${summary_result_text}
        SUMMARY
        exit ${summary_result_code}
    }

    output {
    }
}
