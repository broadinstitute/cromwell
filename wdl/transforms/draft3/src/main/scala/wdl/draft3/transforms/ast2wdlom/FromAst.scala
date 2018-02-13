package wdl.draft3.transforms.ast2wdlom

import wdl.draft3.parser.WdlParser.Ast

trait FromAst[A] extends FromAtoB[Ast, A]
