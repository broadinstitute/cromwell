task summary {
  String bfile
  command {
    ~/plink --bfile ${bfile} --missing --hardy --out foo --allow-no-sex
  }
  output {
    File hwe = "foo.hwe"
    File log = "foo.log"
    File imiss = "foo.imiss"
    File lmiss = "foo.lmiss"
  }
  meta {
    author: "A name"
    email: "aname@broadinstitute.org"
  }
}

workflow test1 {
  File bfile
  call summary {
     input: bfile = bfile
  }
}