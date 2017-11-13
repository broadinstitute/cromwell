import "sub_workflow_decls_import.wdl" as decls

workflow sub_workflow_decls {
  call decls.sub_decls as sudecls { input:
    uninitialized = "initialized",
    wantsOverride = "overridden"
  }

  output {
    String out = sudecls.out
    String init = sudecls.init
    String override = sudecls.override
  }
}
