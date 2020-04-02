#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_bcs.inc.sh" || source test_bcs.inc.sh

cromwell::build::setup_common_environment

cromwell::build::bcs::setup_bcs_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

# Instead of excluding tests, only include a fixed list of tests. This is because due to the
# numerous issues below, contributors did not like having to constantly update the exclude lists.
# https://github.com/broadinstitute/cromwell/issues/3522
# https://github.com/broadinstitute/cromwell/issues/3523
# https://github.com/broadinstitute/cromwell/issues/3524
# https://github.com/broadinstitute/cromwell/issues/3518
# https://github.com/broadinstitute/cromwell/issues/3519
include_tests=( \
    -i abort.instant_abort \
    -i abort.sub_workflow_abort \
    -i aliased_subworkflows \
    -i array_io \
    -i array_literal_locations \
    -i arrays_scatters_ifs \
    -i bad_docker_name \
    -i bad_workflow_failure_mode \
    -i cacheBetweenWF \
    -i cacheWithinWF \
    -i chainfail \
    -i complex_types_files \
    -i composedenginefunctions \
    -i cwl_input_binding_expression \
    -i cwl_optionals \
    -i declarations \
    -i declarations_as_nodes \
    -i declarations_in_ifs \
    -i default_runtime_attributes \
    -i defined_function \
    -i dont_strip_line_prefix \
    -i dot_dir_stuck_running \
    -i draft3_declaration_chain \
    -i draft3_default_input_overrides \
    -i draft3_empty \
    -i draft3_import_structs \
    -i draft3_lots_of_nesting \
    -i draft3_nested_scatter \
    -i draft3_nested_struct \
    -i draft3_passthrough_value \
    -i draft3_sizeenginefunction \
    -i draft3_struct_output \
    -i draft3_taskless_engine_functions \
    -i empty_scatter \
    -i empty_string \
    -i exit \
    -i expression_lib_cwl \
    -i failures.terminal_status \
    -i filearrayoutput \
    -i floating_tags \
    -i forkjoin \
    -i hello_cwl \
    -i if_then_else_expressions \
    -i ifs_upstream_and_downstream \
    -i input_mirror \
    -i invalid_inputs_json \
    -i invalid_labels \
    -i invalid_options_json \
    -i invalid_runtime_attributes \
    -i invalid_wdl \
    -i length \
    # -i long_cmd \ # 2019-08-05 consistently timing out trying to read a < 100KB file in 60 seconds
    -i lots_of_nesting \
    -i member_access \
    -i missing_imports \
    -i missing_sub_inputs \
    -i multiline_command_line \
    -i multiplesourcedarray \
    -i nested_lookups \
    -i null_input_values \
    -i object_access \
    -i optional_declarations \
    -i optional_parameter \
    -i output_filename_interpolation \
    -i output_redirection \
    -i passingfiles \
    -i prefix \
    -i public_http_import \
    -i readFromCacheFalse \
    -i read_tsv \
    -i read_write_json \
    -i read_write_map \
    -i referencingpreviousinputsandoutputs \
    -i runtime_attribute_expressions \
    -i runtime_failOnStderr \
    -i scatterchain \
    -i scattergather \
    -i scatters_in_ifs \
    -i select_functions \
    -i simple_if \
    -i simple_if_workflow_outputs \
    -i single_to_array_coercion \
    -i sizeenginefunction \
    -i square \
    -i stdout_stderr_passing \
    -i string_interpolation \
    -i sub_function \
    -i sub_workflow_decls \
    -i sub_workflow_hello_world \
    -i sub_workflow_interactions \
    -i sub_workflow_interactions_scatter \
    -i sub_workflow_no_output \
    -i sub_workflow_var_refs \
    -i subdirectory \
    -i subworkflows_in_ifs \
    -i taskless_engine_functions \
    -i test_file_outputs_from_input \
    -i three_step__subwf_cwl \
    -i unexpected_call_input_failure \
    -i unexpected_subworkflow_call_input_failure \
    -i valid_labels \
    -i variable_scoping \
    -i wdl_function_locations \
    -i workflow_output_declarations \
    -i workflow_type_and_version_default \
    -i workflow_type_and_version_wdl \
    -i workflowenginefunctions \
    -i writeToCache \
    -i write_lines \
    -i write_lines_files \
    -i write_tsv \
)

cromwell::build::run_centaur \
    -p 100 \
    -t 1m \
    -e localdockertest \
    "${include_tests[@]}" \

cromwell::build::generate_code_coverage
