
package cromwell.parser;
import java.util.*;
import java.io.IOException;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import org.apache.commons.codec.binary.Base64;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.lang.reflect.Method;
public class WdlParser {
    private static Map<Integer, List<TerminalIdentifier>> nonterminal_first;
    private static Map<Integer, List<TerminalIdentifier>> nonterminal_follow;
    private static Map<Integer, List<TerminalIdentifier>> rule_first;
    private static Map<Integer, List<String>> nonterminal_rules;
    private static Map<Integer, String> rules;
    public static WdlTerminalMap terminal_map = new WdlTerminalMap(WdlTerminalIdentifier.values());
    public WdlParser() {
        try {
            lexer_init();
        } catch(Exception e) {}
    }
    public static String join(Collection<?> s, String delimiter) {
        StringBuilder builder = new StringBuilder();
        Iterator iter = s.iterator();
        while (iter.hasNext()) {
            builder.append(iter.next());
            if (!iter.hasNext()) {
                break;
            }
            builder.append(delimiter);
        }
        return builder.toString();
    }
    public static String getIndentString(int spaces) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < spaces; i++) {
            sb.append(' ');
        }
        return sb.toString();
    }
    public static String readStdin() throws IOException {
        InputStreamReader stream = new InputStreamReader(System.in, "utf-8");
        char buffer[] = new char[System.in.available()];
        try {
            stream.read(buffer, 0, System.in.available());
        } finally {
            stream.close();
        }
        return new String(buffer);
    }
    public static String readFile(String path) throws IOException {
        FileInputStream stream = new FileInputStream(new File(path));
        try {
            FileChannel fc = stream.getChannel();
            MappedByteBuffer bb = fc.map(FileChannel.MapMode.READ_ONLY, 0, fc.size());
            /* Instead of using default, pass in a decoder. */
            return Charset.defaultCharset().decode(bb).toString();
        }
        finally {
            stream.close();
        }
    }
    public static class SyntaxError extends Exception {
        public SyntaxError(String message) {
            super(message);
        }
    }
    public interface SyntaxErrorFormatter {
        /* Called when the parser runs out of tokens but isn't finished parsing. */
        String unexpectedEof(String method, List<TerminalIdentifier> expected, List<String> nt_rules);
        /* Called when the parser finished parsing but there are still tokens left in the stream. */
        String excessTokens(String method, Terminal terminal);
        /* Called when the parser is expecting one token and gets another. */
        String unexpectedSymbol(String method, Terminal actual, List<TerminalIdentifier> expected, String rule);
        /* Called when the parser is expecing a tokens but there are no more tokens. */
        String noMoreTokens(String method, TerminalIdentifier expecting, Terminal last);
        /* Invalid terminal is found in the token stream. */
        String invalidTerminal(String method, Terminal invalid);
    }
    public static class TokenStream extends ArrayList<Terminal> {
        private int index;
        public TokenStream(List<Terminal> terminals) {
            super(terminals);
            reset();
        }
        public TokenStream() {
            reset();
        }
        public void reset() {
            this.index = 0;
        }
        public Terminal advance() {
            this.index += 1;
            return this.current();
        }
        public Terminal current() {
            try {
                return this.get(this.index);
            } catch (IndexOutOfBoundsException e) {
                return null;
            }
        }
        public Terminal last() {
          return this.get(this.size() - 1);
        }
    }
    public static class NonTerminal {
        private int id;
        private String string;
        NonTerminal(int id, String string) {
            this.id = id;
            this.string = string;
        }
        public int getId() {
            return this.id;
        }
        public String getString() {
            return this.string;
        }
        public String toString() {
            return this.string;
        }
    }
    public interface AstTransform {}
    public static class AstTransformNodeCreator implements AstTransform {
        private String name;
        private LinkedHashMap<String, Integer> parameters;
        AstTransformNodeCreator(String name, LinkedHashMap<String, Integer> parameters) {
            this.name = name;
            this.parameters = parameters;
        }
        public Map<String, Integer> getParameters() {
            return this.parameters;
        }
        public String getName() {
            return this.name;
        }
        public String toString() {
            LinkedList<String> items = new LinkedList<String>();
            for (final Map.Entry<String, Integer> entry : this.parameters.entrySet()) {
                items.add(entry.getKey() + "=$" + entry.getValue().toString());
            }
            return "AstNodeCreator: " + this.name + "( " + join(items, ", ") + " )";
        }
    }
    public static class AstTransformSubstitution implements AstTransform {
        private int index;
        AstTransformSubstitution(int index) {
            this.index = index;
        }
        public int getIndex() {
            return this.index;
        }
        public String toString() {
            return "AstSubstitution: $" + Integer.toString(this.index);
        }
    }
    public interface AstNode {
        public String toString();
        public String toPrettyString();
        public String toPrettyString(int indent);
    }
    public static class AstList extends ArrayList<AstNode> implements AstNode {
        public String toString() {
            return "[" + join(this, ", ") + "]";
        }
        public String toPrettyString() {
            return toPrettyString(0);
        }
        public String toPrettyString(int indent) {
            String spaces = getIndentString(indent);
            if (this.size() == 0) {
                return spaces + "[]";
            }
            ArrayList<String> elements = new ArrayList<String>();
            for ( AstNode node : this ) {
                elements.add(node.toPrettyString(indent + 2));
            }
            return spaces + "[\n" + join(elements, ",\n") + "\n" + spaces + "]";
        }
    }
    public static class Ast implements AstNode {
        private String name;
        private Map<String, AstNode> attributes;
        Ast(String name, Map<String, AstNode> attributes) {
            this.name = name;
            this.attributes = attributes;
        }
        public AstNode getAttribute(String name) {
            return this.attributes.get(name);
        }
        public Map<String, AstNode> getAttributes() {
            return this.attributes;
        }
        public String getName() {
            return this.name;
        }
        public String toString() {
            Formatter formatter = new Formatter(new StringBuilder(), Locale.US);
            LinkedList<String> attributes = new LinkedList<String>();
            for (final Map.Entry<String, AstNode> attribute : this.attributes.entrySet()) {
                final String name = attribute.getKey();
                final AstNode node = attribute.getValue();
                final String nodeStr = (node == null) ? "None" : node.toString();
                attributes.add(name + "=" + nodeStr);
            }
            formatter.format("(%s: %s)", this.name, join(attributes, ", "));
            return formatter.toString();
        }
        public String toPrettyString() {
            return toPrettyString(0);
        }
        public String toPrettyString(int indent) {
            String spaces = getIndentString(indent);
            ArrayList<String> children = new ArrayList<String>();
            for( Map.Entry<String, AstNode> attribute : this.attributes.entrySet() ) {
                String valueString = attribute.getValue() == null ? "None" : attribute.getValue().toPrettyString(indent + 2).trim();
                children.add(spaces + "  " + attribute.getKey() + "=" + valueString);
            }
            return spaces + "(" + this.name + ":\n" + join(children, ",\n") + "\n" + spaces + ")";
        }
    }
    public interface ParseTreeNode {
        public AstNode toAst();
        public String toString();
        public String toPrettyString();
        public String toPrettyString(int indent);
    }
    public static class Terminal implements AstNode, ParseTreeNode
    {
        private int id;
        private String terminal_str;
        private String source_string;
        private String resource;
        private int line;
        private int col;
        public Terminal(int id, String terminal_str, String source_string, String resource, int line, int col) {
            this.id = id;
            this.terminal_str = terminal_str;
            this.source_string = source_string;
            this.resource = resource;
            this.line = line;
            this.col = col;
        }
        public int getId() {
            return this.id;
        }
        public String getTerminalStr() {
            return this.terminal_str;
        }
        public String getSourceString() {
            return this.source_string;
        }
        public String getResource() {
            return this.resource;
        }
        public int getLine() {
            return this.line;
        }
        public int getColumn() {
            return this.col;
        }
        public String toString() {
            byte[] source_string_bytes;
            try {
                source_string_bytes = this.getSourceString().getBytes("UTF-8");
            } catch (java.io.UnsupportedEncodingException e) {
                source_string_bytes = this.getSourceString().getBytes();
            }
            return String.format("<%s:%d:%d %s \"%s\">",
                this.getResource(),
                this.getLine(),
                this.getColumn(),
                this.getTerminalStr(),
                Base64.encodeBase64String(source_string_bytes)
            );
        }
        public String toPrettyString() {
            return toPrettyString(0);
        }
        public String toPrettyString(int indent) {
            return getIndentString(indent) + this.toString();
        }
        public AstNode toAst() { return this; }
    }
    public static class ParseTree implements ParseTreeNode {
        private NonTerminal nonterminal;
        private ArrayList<ParseTreeNode> children;
        private boolean isExpr, isNud, isPrefix, isInfix, isExprNud;
        private int nudMorphemeCount;
        private Terminal listSeparator;
        private boolean list;
        private AstTransform astTransform;
        ParseTree(NonTerminal nonterminal) {
            this.nonterminal = nonterminal;
            this.children = new ArrayList<ParseTreeNode>();
            this.astTransform = null;
            this.isExpr = false;
            this.isNud = false;
            this.isPrefix = false;
            this.isInfix = false;
            this.isExprNud = false;
            this.nudMorphemeCount = 0;
            this.listSeparator = null;
            this.list = false;
        }
        public void setExpr(boolean value) { this.isExpr = value; }
        public void setNud(boolean value) { this.isNud = value; }
        public void setPrefix(boolean value) { this.isPrefix = value; }
        public void setInfix(boolean value) { this.isInfix = value; }
        public void setExprNud(boolean value) { this.isExprNud = value; }
        public void setAstTransformation(AstTransform value) { this.astTransform = value; }
        public void setNudMorphemeCount(int value) { this.nudMorphemeCount = value; }
        public void setList(boolean value) { this.list = value; }
        public void setListSeparator(Terminal value) { this.listSeparator = value; }
        public int getNudMorphemeCount() { return this.nudMorphemeCount; }
        public List<ParseTreeNode> getChildren() { return this.children; }
        public boolean isInfix() { return this.isInfix; }
        public boolean isPrefix() { return this.isPrefix; }
        public boolean isExpr() { return this.isExpr; }
        public boolean isNud() { return this.isNud; }
        public boolean isExprNud() { return this.isExprNud; }
        public void add(ParseTreeNode tree) {
            if (this.children == null) {
                this.children = new ArrayList<ParseTreeNode>();
            }
            this.children.add(tree);
        }
        private boolean isCompoundNud() {
            if ( this.children.size() > 0 && this.children.get(0) instanceof ParseTree ) {
                ParseTree child = (ParseTree) this.children.get(0);
                if ( child.isNud() && !child.isPrefix() && !this.isExprNud() && !this.isInfix() ) {
                    return true;
                }
            }
            return false;
        }
        public AstNode toAst() {
            if ( this.list == true ) {
                AstList astList = new AstList();
                int end = this.children.size() - 1;
                if ( this.children.size() == 0 ) {
                    return astList;
                }
                for (int i = 0; i < this.children.size() - 1; i++) {
                    if (this.children.get(i) instanceof Terminal && this.listSeparator != null && ((Terminal)this.children.get(i)).id == this.listSeparator.id)
                        continue;
                    astList.add(this.children.get(i).toAst());
                }
                astList.addAll((AstList) this.children.get(this.children.size() - 1).toAst());
                return astList;
            } else if ( this.isExpr ) {
                if ( this.astTransform instanceof AstTransformSubstitution ) {
                    AstTransformSubstitution astSubstitution = (AstTransformSubstitution) astTransform;
                    return this.children.get(astSubstitution.getIndex()).toAst();
                } else if ( this.astTransform instanceof AstTransformNodeCreator ) {
                    AstTransformNodeCreator astNodeCreator = (AstTransformNodeCreator) this.astTransform;
                    LinkedHashMap<String, AstNode> parameters = new LinkedHashMap<String, AstNode>();
                    ParseTreeNode child;
                    for ( final Map.Entry<String, Integer> parameter : astNodeCreator.getParameters().entrySet() ) {
                        String name = parameter.getKey();
                        int index = parameter.getValue().intValue();
                        if ( index == '$' ) {
                            child = this.children.get(0);
                        } else if ( this.isCompoundNud() ) {
                            ParseTree firstChild = (ParseTree) this.children.get(0);
                            if ( index < firstChild.getNudMorphemeCount() ) {
                                child = firstChild.getChildren().get(index);
                            } else {
                                index = index - firstChild.getNudMorphemeCount() + 1;
                                child = this.children.get(index);
                            }
                        } else if ( this.children.size() == 1 && !(this.children.get(0) instanceof ParseTree) && !(this.children.get(0) instanceof List) ) {
                            // TODO: I don't think this should ever be called
                            child = this.children.get(0);
                        } else {
                            child = this.children.get(index);
                        }
                        parameters.put(name, child.toAst());
                    }
                    return new Ast(astNodeCreator.getName(), parameters);
                }
            } else {
                AstTransformSubstitution defaultAction = new AstTransformSubstitution(0);
                AstTransform action = this.astTransform != null ? this.astTransform : defaultAction;
                if (this.children.size() == 0) return null;
                if (action instanceof AstTransformSubstitution) {
                    AstTransformSubstitution astSubstitution = (AstTransformSubstitution) action;
                    return this.children.get(astSubstitution.getIndex()).toAst();
                } else if (action instanceof AstTransformNodeCreator) {
                    AstTransformNodeCreator astNodeCreator = (AstTransformNodeCreator) action;
                    LinkedHashMap<String, AstNode> evaluatedParameters = new LinkedHashMap<String, AstNode>();
                    for ( Map.Entry<String, Integer> baseParameter : astNodeCreator.getParameters().entrySet() ) {
                        String name = baseParameter.getKey();
                        int index2 = baseParameter.getValue().intValue();
                        evaluatedParameters.put(name, this.children.get(index2).toAst());
                    }
                    return new Ast(astNodeCreator.getName(), evaluatedParameters);
                }
            }
            return null;
        }
        public String toString() {
          ArrayList<String> children = new ArrayList<String>();
          for (ParseTreeNode child : this.children) {
            children.add(child.toString());
          }
          return "(" + this.nonterminal.getString() + ": " + join(children, ", ") + ")";
        }
        public String toPrettyString() {
          return toPrettyString(0);
        }
        public String toPrettyString(int indent) {
          if (this.children.size() == 0) {
            return "(" + this.nonterminal.toString() + ": )";
          }
          String spaces = getIndentString(indent);
          ArrayList<String> children = new ArrayList<String>();
          for ( ParseTreeNode node : this.children ) {
            String sub = node.toPrettyString(indent + 2).trim();
            children.add(spaces + "  " +  sub);
          }
          return spaces + "(" + this.nonterminal.toString() + ":\n" + join(children, ",\n") + "\n" + spaces + ")";
        }
    }
    private static class ParserContext {
        public TokenStream tokens;
        public SyntaxErrorFormatter error_formatter;
        public String nonterminal;
        public String rule;
        public ParserContext(TokenStream tokens, SyntaxErrorFormatter error_formatter) {
            this.tokens = tokens;
            this.error_formatter = error_formatter;
        }
    }
    private static class DefaultSyntaxErrorFormatter implements SyntaxErrorFormatter {
        public String unexpectedEof(String method, List<TerminalIdentifier> expected, List<String> nt_rules) {
            return "Error: unexpected end of file";
        }
        public String excessTokens(String method, Terminal terminal) {
            return "Finished parsing without consuming all tokens.";
        }
        public String unexpectedSymbol(String method, Terminal actual, List<TerminalIdentifier> expected, String rule) {
            ArrayList<String> expected_terminals = new ArrayList<String>();
            for ( TerminalIdentifier e : expected ) {
                expected_terminals.add(e.string());
            }
            return String.format(
                "Unexpected symbol (line %d, col %d) when parsing parse_%s.  Expected %s, got %s.",
                actual.getLine(), actual.getColumn(), method, join(expected_terminals, ", "), actual.toPrettyString()
            );
        }
        public String noMoreTokens(String method, TerminalIdentifier expecting, Terminal last) {
            return "No more tokens.  Expecting " + expecting.string();
        }
        public String invalidTerminal(String method, Terminal invalid) {
            return "Invalid symbol ID: "+invalid.getId()+" ("+invalid.getTerminalStr()+")";
        }
    }
    public interface TerminalMap {
        TerminalIdentifier get(String string);
        TerminalIdentifier get(int id);
        boolean isValid(String string);
        boolean isValid(int id);
    }
    public static class WdlTerminalMap implements TerminalMap {
        private Map<Integer, TerminalIdentifier> id_to_term;
        private Map<String, TerminalIdentifier> str_to_term;
        WdlTerminalMap(WdlTerminalIdentifier[] terminals) {
            id_to_term = new HashMap<Integer, TerminalIdentifier>();
            str_to_term = new HashMap<String, TerminalIdentifier>();
            for( WdlTerminalIdentifier terminal : terminals ) {
                Integer id = new Integer(terminal.id());
                String str = terminal.string();
                id_to_term.put(id, terminal);
                str_to_term.put(str, terminal);
            }
        }
        public TerminalIdentifier get(String string) { return this.str_to_term.get(string); }
        public TerminalIdentifier get(int id) { return this.id_to_term.get(id); }
        public boolean isValid(String string) { return this.str_to_term.containsKey(string); }
        public boolean isValid(int id) { return this.id_to_term.containsKey(id); }
    }
    public interface TerminalIdentifier {
        public int id();
        public String string();
    }
    public enum WdlTerminalIdentifier implements TerminalIdentifier {
        TERMINAL_RAW_CMD_END(0, "raw_cmd_end"),
        TERMINAL_WORKFLOW(1, "workflow"),
        TERMINAL_INPUT(2, "input"),
        TERMINAL_IDENTIFIER(3, "identifier"),
        TERMINAL_LBRACE(4, "lbrace"),
        TERMINAL_DOT(5, "dot"),
        TERMINAL_COMMA(6, "comma"),
        TERMINAL_CMD_PARAM_START(7, "cmd_param_start"),
        TERMINAL_ASTERISK(8, "asterisk"),
        TERMINAL_GT(9, "gt"),
        TERMINAL_RAW_COMMAND(10, "raw_command"),
        TERMINAL_RPAREN(11, "rparen"),
        TERMINAL_STRING(12, "string"),
        TERMINAL_RAW_CMD_START(13, "raw_cmd_start"),
        TERMINAL_PARAMETER_META(14, "parameter_meta"),
        TERMINAL_LT(15, "lt"),
        TERMINAL_NOT(16, "not"),
        TERMINAL_SQUOTE_STRING(17, "squote_string"),
        TERMINAL_QMARK(18, "qmark"),
        TERMINAL_IN(19, "in"),
        TERMINAL_LPAREN(20, "lparen"),
        TERMINAL_NS_IDENTIFIER(21, "ns_identifier"),
        TERMINAL_SLASH(22, "slash"),
        TERMINAL_DOUBLE_EQUAL(23, "double_equal"),
        TERMINAL_CMD_PART(24, "cmd_part"),
        TERMINAL_OUTPUT(25, "output"),
        TERMINAL_AS(26, "as"),
        TERMINAL_META(27, "meta"),
        TERMINAL_LSQUARE(28, "lsquare"),
        TERMINAL_PLUS(29, "plus"),
        TERMINAL_RSQUARE(30, "rsquare"),
        TERMINAL_OBJECT(31, "object"),
        TERMINAL_DQUOTE_STRING(32, "dquote_string"),
        TERMINAL_NOT_EQUAL(33, "not_equal"),
        TERMINAL_INTEGER(34, "integer"),
        TERMINAL_WHILE(35, "while"),
        TERMINAL_TASK(36, "task"),
        TERMINAL_CALL(37, "call"),
        TERMINAL_DASH(38, "dash"),
        TERMINAL_TYPE(39, "type"),
        TERMINAL_FLOAT(40, "float"),
        TERMINAL_IMPORT(41, "import"),
        TERMINAL_E(42, "e"),
        TERMINAL_DOUBLE_AMPERSAND(43, "double_ampersand"),
        TERMINAL_DOUBLE_PIPE(44, "double_pipe"),
        TERMINAL_IF(45, "if"),
        TERMINAL_SCATTER(46, "scatter"),
        TERMINAL_BOOLEAN(47, "boolean"),
        TERMINAL_RUNTIME(48, "runtime"),
        TERMINAL_RBRACE(49, "rbrace"),
        TERMINAL_COLON(50, "colon"),
        TERMINAL_CMD_PARAM_END(51, "cmd_param_end"),
        TERMINAL_TYPE_E(52, "type_e"),
        TERMINAL_EQUAL(53, "equal"),
        TERMINAL_GTEQ(54, "gteq"),
        TERMINAL_CMD_ATTR_HINT(55, "cmd_attr_hint"),
        TERMINAL_PERCENT(56, "percent"),
        TERMINAL_LTEQ(57, "lteq"),
        END_SENTINAL(-3, "END_SENTINAL");
        private final int id;
        private final String string;
        WdlTerminalIdentifier(int id, String string) {
            this.id = id;
            this.string = string;
        }
        public int id() {return id;}
        public String string() {return string;}
    }
    /* table[nonterminal][terminal] = rule */
    private static final int[][] table = {
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 58, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 40, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 41, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 41, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 43, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 42, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 42, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 19, -1, -1, -1, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 19, -1, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 19, 22, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 61, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, 18, -1, -1, -1, -1, -1, -1, -1, 18, -1, -1, -1, 18, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 18, -1, 18, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, -1, -1, -1, -1, -1, -1, -1, -1, 18, -1, -1, -1, 15, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 99, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 66, -1, 66, -1, 66, -1, -1, -1, -1, -1, 66, 66, -1, -1, 69, -1, -1, 66, -1, -1, -1, -1, -1 },
        { -1, -1, 82, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 83, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, 46, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 80, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 60, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 63, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, 90, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, 89, 86, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 89, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 8, -1, -1, -1, -1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, 131, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 132, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 47, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, 78, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 79, -1, 79, -1, 79, -1, -1, -1, -1, -1, 79, 79, -1, -1, 79, -1, -1, 79, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, 64, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 95, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 24, -1, -1, -1, 27, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 25, -1, 28, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 26, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 52, -1, -1, -1, -1, -1, -1, -1, -1, -1, 53, -1, -1, 52, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 59, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, 88, -1, -1, -1, 87, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 88, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, 102, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 103, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 39, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 36, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 73, -1, 71, -1, 72, -1, -1, -1, -1, -1, 74, 75, -1, -1, -1, -1, -1, 72, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 38, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 37, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 67, -1, 67, -1, 67, -1, -1, -1, -1, -1, 67, 67, -1, -1, 68, -1, -1, 67, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 98, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 98, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 20, -1, -1, -1, 20, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 20, -1, 20, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 20, 21, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, 12, -1, -1, -1, -1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, 77, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 76, -1, -1, -1, -1, -1, -1, -1, -1, 77, -1, 77, -1, 77, -1, -1, -1, -1, -1, 77, 77, -1, -1, 77, -1, -1, 77, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, 35, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 34, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 16, -1, -1, -1, -1, -1, -1, -1, -1, 17, -1, -1, -1, 16, -1, -1, -1, -1, -1 },
        { -1, -1, 97, -1, -1, -1, -1, -1, -1, -1, 97, -1, -1, -1, 97, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 97, -1, 97, -1, -1, -1, -1, -1, -1, -1, 97, -1, 97, -1, 97, -1, -1, -1, -1, -1, 97, 97, -1, 97, 97, -1, -1, 97, 96, -1, -1, -1, -1 },
        { -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 57, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 104, -1, -1, -1, -1, -1, -1, -1, -1, 101, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 101, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 130, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 133, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 32, -1, -1, -1, -1, -1, -1, 29, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 29, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, 50, -1, -1, -1, -1, -1, -1, -1, -1, -1, 48, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 49, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, 81, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 84, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, 124, -1, -1, -1, -1, 125, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 125, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 55, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 70, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 14, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 56, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 56, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, 85, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 91, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 92, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 65, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 33, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 93, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 51, -1, -1, -1, -1, -1, -1, -1, -1, -1, 54, -1, -1, 51, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 94, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 100, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, 123, -1, -1, -1, -1, -1, -1, -1, 126, 123, -1, -1, -1, 123, 123, -1, -1, 123, -1, -1, -1, -1, -1, -1, -1, 123, 123, 126, 123, 123, -1, 123, -1, -1, -1, 123, -1, 123, -1, 123, -1, -1, -1, -1, 123, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { 31, -1, -1, -1, -1, -1, -1, 30, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 30, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1, 44, -1, -1, -1, -1, -1, -1, -1, -1, -1, 44, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 44, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 45, -1, -1, -1, -1, -1, -1 },
    };
    static {
        Map<Integer, List<TerminalIdentifier>> map = new HashMap<Integer, List<TerminalIdentifier>>();
        map.put(58, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(59, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(60, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(61, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(62, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(63, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(64, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_EQUAL,
        }));
        map.put(65, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(66, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
        }));
        map.put(67, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
        }));
        map.put(68, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(69, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
        }));
        map.put(70, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(71, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
        }));
        map.put(72, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
        }));
        map.put(73, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(74, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(75, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(76, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
        }));
        map.put(77, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(78, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LBRACE,
        }));
        map.put(79, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LBRACE,
        }));
        map.put(80, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_SCATTER,
        }));
        map.put(81, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(82, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(83, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
        }));
        map.put(84, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(85, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_META,
        }));
        map.put(86, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(87, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(88, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
        }));
        map.put(89, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
        }));
        map.put(90, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
        }));
        map.put(91, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(92, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(93, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(94, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(95, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(96, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(97, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(98, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(99, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_EQUAL,
        }));
        map.put(100, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(101, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
        }));
        map.put(102, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(103, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(104, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(105, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_ASTERISK,
            WdlTerminalIdentifier.TERMINAL_QMARK,
            WdlTerminalIdentifier.TERMINAL_PLUS,
        }));
        map.put(106, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
        }));
        map.put(107, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(108, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
        }));
        map.put(109, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
        }));
        map.put(110, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(111, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(112, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(113, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LBRACE,
        }));
        map.put(114, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(115, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(116, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(117, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
        }));
        map.put(118, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
        }));
        map.put(119, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(120, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(121, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IF,
        }));
        map.put(122, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(123, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_STRING,
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
        }));
        map.put(124, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(125, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_ASTERISK,
            WdlTerminalIdentifier.TERMINAL_QMARK,
            WdlTerminalIdentifier.TERMINAL_PLUS,
        }));
        nonterminal_first = Collections.unmodifiableMap(map);
    }
    static {
        Map<Integer, List<TerminalIdentifier>> map = new HashMap<Integer, List<TerminalIdentifier>>();
        map.put(58, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(59, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(60, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(61, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(62, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(63, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(64, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_META,
        }));
        map.put(65, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(66, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(67, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_RAW_CMD_END,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(68, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(69, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_IMPORT,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(70, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(71, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(72, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(73, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(74, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(75, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(76, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(77, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_RSQUARE,
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_NOT_EQUAL,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_COMMA,
            WdlTerminalIdentifier.TERMINAL_GT,
            WdlTerminalIdentifier.TERMINAL_ASTERISK,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_RPAREN,
            WdlTerminalIdentifier.TERMINAL_STRING,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
            WdlTerminalIdentifier.TERMINAL_LT,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_DOUBLE_AMPERSAND,
            WdlTerminalIdentifier.TERMINAL_DOUBLE_PIPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SLASH,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_DOUBLE_EQUAL,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_GTEQ,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
            WdlTerminalIdentifier.TERMINAL_PERCENT,
            WdlTerminalIdentifier.TERMINAL_LTEQ,
        }));
        map.put(78, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(79, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(80, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(81, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(82, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_RSQUARE,
        }));
        map.put(83, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(84, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(85, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(86, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(87, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RSQUARE,
        }));
        map.put(88, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(89, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(90, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(91, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(92, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(93, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(94, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_IMPORT,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(95, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(96, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_LBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(97, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_RAW_CMD_END,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(98, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(99, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_META,
        }));
        map.put(100, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(101, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(102, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RSQUARE,
        }));
        map.put(103, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(104, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_CMD_END,
        }));
        map.put(105, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_END,
        }));
        map.put(106, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(107, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RPAREN,
            WdlTerminalIdentifier.TERMINAL_RSQUARE,
        }));
        map.put(108, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(109, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(110, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(111, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_IMPORT,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(112, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(113, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(114, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
            WdlTerminalIdentifier.TERMINAL_INPUT,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(115, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_LBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(116, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(117, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(118, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(119, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(120, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(121, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(122, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
            WdlTerminalIdentifier.TERMINAL_RBRACE,
        }));
        map.put(123, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RPAREN,
            WdlTerminalIdentifier.TERMINAL_RSQUARE,
        }));
        map.put(124, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_CMD_END,
        }));
        map.put(125, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_END,
        }));
        nonterminal_follow = Collections.unmodifiableMap(map);
    }
    static {
        Map<Integer, List<TerminalIdentifier>> map = new HashMap<Integer, List<TerminalIdentifier>>();
        map.put(0, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
        }));
        map.put(1, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
        }));
        map.put(2, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(3, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(4, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(5, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(6, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(7, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(8, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
            WdlTerminalIdentifier.TERMINAL_IMPORT,
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(9, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(10, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
        }));
        map.put(11, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(12, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(13, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IMPORT,
        }));
        map.put(14, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(15, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(16, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(17, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(18, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(19, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(20, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
            WdlTerminalIdentifier.TERMINAL_META,
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(21, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(22, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(23, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TASK,
        }));
        map.put(24, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
        }));
        map.put(25, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
        }));
        map.put(26, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
        }));
        map.put(27, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(28, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_META,
        }));
        map.put(29, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(30, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(31, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(32, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(33, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
        }));
        map.put(34, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PART,
        }));
        map.put(35, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
        }));
        map.put(36, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
        }));
        map.put(37, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
        }));
        map.put(38, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(39, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(40, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(41, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(42, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(43, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(44, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_ASTERISK,
            WdlTerminalIdentifier.TERMINAL_QMARK,
            WdlTerminalIdentifier.TERMINAL_PLUS,
        }));
        map.put(45, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(46, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
        }));
        map.put(47, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
        }));
        map.put(48, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_QMARK,
        }));
        map.put(49, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PLUS,
        }));
        map.put(50, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_ASTERISK,
        }));
        map.put(51, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(52, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(53, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(54, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(55, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OUTPUT,
        }));
        map.put(56, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(57, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_RUNTIME,
        }));
        map.put(58, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
        }));
        map.put(59, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_META,
        }));
        map.put(60, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(61, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(62, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(63, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(64, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LBRACE,
        }));
        map.put(65, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(66, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
        }));
        map.put(67, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
            WdlTerminalIdentifier.TERMINAL_CALL,
            WdlTerminalIdentifier.TERMINAL_IF,
            WdlTerminalIdentifier.TERMINAL_SCATTER,
        }));
        map.put(68, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(69, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(70, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WORKFLOW,
        }));
        map.put(71, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(72, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(73, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
        }));
        map.put(74, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IF,
        }));
        map.put(75, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_SCATTER,
        }));
        map.put(76, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(77, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(78, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LBRACE,
        }));
        map.put(79, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(80, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_CALL,
        }));
        map.put(81, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
        }));
        map.put(82, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
        }));
        map.put(83, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(84, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(85, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LBRACE,
        }));
        map.put(86, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(87, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(88, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(89, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(90, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INPUT,
        }));
        map.put(91, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(92, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_AS,
        }));
        map.put(93, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_WHILE,
        }));
        map.put(94, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IF,
        }));
        map.put(95, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_SCATTER,
        }));
        map.put(96, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_EQUAL,
        }));
        map.put(97, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(98, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(99, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_EQUAL,
        }));
        map.put(100, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(101, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE_E,
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(102, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(103, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(104, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(105, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(106, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_TYPE,
        }));
        map.put(107, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(108, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(109, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(110, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(111, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(112, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(113, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(114, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(115, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(116, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(117, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(118, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(119, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(120, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
        }));
        map.put(121, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_PLUS,
        }));
        map.put(122, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_DASH,
        }));
        map.put(123, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_NOT,
            WdlTerminalIdentifier.TERMINAL_FLOAT,
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
            WdlTerminalIdentifier.TERMINAL_E,
            WdlTerminalIdentifier.TERMINAL_PLUS,
            WdlTerminalIdentifier.TERMINAL_LPAREN,
            WdlTerminalIdentifier.TERMINAL_OBJECT,
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
            WdlTerminalIdentifier.TERMINAL_INTEGER,
            WdlTerminalIdentifier.TERMINAL_DASH,
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(124, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(125, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(126, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(127, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(128, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(129, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(130, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(131, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_COMMA,
        }));
        map.put(132, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(133, Arrays.asList(new TerminalIdentifier[] {
        }));
        map.put(134, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_OBJECT,
        }));
        map.put(135, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LSQUARE,
        }));
        map.put(136, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_LPAREN,
        }));
        map.put(137, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_STRING,
        }));
        map.put(138, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
        }));
        map.put(139, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_BOOLEAN,
        }));
        map.put(140, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_INTEGER,
        }));
        map.put(141, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_FLOAT,
        }));
        map.put(142, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING,
        }));
        map.put(143, Arrays.asList(new TerminalIdentifier[] {
            WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING,
        }));
        rule_first = Collections.unmodifiableMap(map);
    }
    static {
        Map<Integer, List<String>> map = new HashMap<Integer, List<String>>();
        map.put(58, new ArrayList<String>());
        map.put(59, new ArrayList<String>());
        map.put(60, new ArrayList<String>());
        map.put(61, new ArrayList<String>());
        map.put(62, new ArrayList<String>());
        map.put(63, new ArrayList<String>());
        map.put(64, new ArrayList<String>());
        map.put(65, new ArrayList<String>());
        map.put(66, new ArrayList<String>());
        map.put(67, new ArrayList<String>());
        map.put(68, new ArrayList<String>());
        map.put(69, new ArrayList<String>());
        map.put(70, new ArrayList<String>());
        map.put(71, new ArrayList<String>());
        map.put(72, new ArrayList<String>());
        map.put(73, new ArrayList<String>());
        map.put(74, new ArrayList<String>());
        map.put(75, new ArrayList<String>());
        map.put(76, new ArrayList<String>());
        map.put(77, new ArrayList<String>());
        map.put(78, new ArrayList<String>());
        map.put(79, new ArrayList<String>());
        map.put(80, new ArrayList<String>());
        map.put(81, new ArrayList<String>());
        map.put(82, new ArrayList<String>());
        map.put(83, new ArrayList<String>());
        map.put(84, new ArrayList<String>());
        map.put(85, new ArrayList<String>());
        map.put(86, new ArrayList<String>());
        map.put(87, new ArrayList<String>());
        map.put(88, new ArrayList<String>());
        map.put(89, new ArrayList<String>());
        map.put(90, new ArrayList<String>());
        map.put(91, new ArrayList<String>());
        map.put(92, new ArrayList<String>());
        map.put(93, new ArrayList<String>());
        map.put(94, new ArrayList<String>());
        map.put(95, new ArrayList<String>());
        map.put(96, new ArrayList<String>());
        map.put(97, new ArrayList<String>());
        map.put(98, new ArrayList<String>());
        map.put(99, new ArrayList<String>());
        map.put(100, new ArrayList<String>());
        map.put(101, new ArrayList<String>());
        map.put(102, new ArrayList<String>());
        map.put(103, new ArrayList<String>());
        map.put(104, new ArrayList<String>());
        map.put(105, new ArrayList<String>());
        map.put(106, new ArrayList<String>());
        map.put(107, new ArrayList<String>());
        map.put(108, new ArrayList<String>());
        map.put(109, new ArrayList<String>());
        map.put(110, new ArrayList<String>());
        map.put(111, new ArrayList<String>());
        map.put(112, new ArrayList<String>());
        map.put(113, new ArrayList<String>());
        map.put(114, new ArrayList<String>());
        map.put(115, new ArrayList<String>());
        map.put(116, new ArrayList<String>());
        map.put(117, new ArrayList<String>());
        map.put(118, new ArrayList<String>());
        map.put(119, new ArrayList<String>());
        map.put(120, new ArrayList<String>());
        map.put(121, new ArrayList<String>());
        map.put(122, new ArrayList<String>());
        map.put(123, new ArrayList<String>());
        map.put(124, new ArrayList<String>());
        map.put(125, new ArrayList<String>());
        map.get(83).add("$_gen0 = $import $_gen1");
        map.get(71).add("$_gen1 = $import $_gen1");
        map.get(71).add("$_gen1 = :_empty");
        map.get(83).add("$_gen0 = :_empty");
        map.get(95).add("$_gen2 = $workflow_or_task $_gen3");
        map.get(100).add("$_gen3 = $workflow_or_task $_gen3");
        map.get(100).add("$_gen3 = :_empty");
        map.get(95).add("$_gen2 = :_empty");
        map.get(74).add("$document = $_gen0 $_gen2 -> Document( imports=$0, definitions=$1 )");
        map.get(119).add("$workflow_or_task = $workflow");
        map.get(119).add("$workflow_or_task = $task");
        map.get(94).add("$_gen4 = $import_namespace");
        map.get(94).add("$_gen4 = :_empty");
        map.get(69).add("$import = :import :string $_gen4 -> Import( uri=$1, namespace=$2 )");
        map.get(111).add("$import_namespace = :as :identifier -> $1");
        map.get(63).add("$_gen5 = $declaration $_gen6");
        map.get(98).add("$_gen6 = $declaration $_gen6");
        map.get(98).add("$_gen6 = :_empty");
        map.get(63).add("$_gen5 = :_empty");
        map.get(61).add("$_gen7 = $sections $_gen8");
        map.get(93).add("$_gen8 = $sections $_gen8");
        map.get(93).add("$_gen8 = :_empty");
        map.get(61).add("$_gen7 = :_empty");
        map.get(108).add("$task = :task :identifier :lbrace $_gen5 $_gen7 :rbrace -> Task( name=$1, declarations=$3, sections=$4 )");
        map.get(81).add("$sections = $command");
        map.get(81).add("$sections = $outputs");
        map.get(81).add("$sections = $runtime");
        map.get(81).add("$sections = $parameter_meta");
        map.get(81).add("$sections = $meta");
        map.get(104).add("$_gen9 = $command_part $_gen10");
        map.get(124).add("$_gen10 = $command_part $_gen10");
        map.get(124).add("$_gen10 = :_empty");
        map.get(104).add("$_gen9 = :_empty");
        map.get(117).add("$command = :raw_command :raw_cmd_start $_gen9 :raw_cmd_end -> RawCommand( parts=$2 )");
        map.get(97).add("$command_part = :cmd_part");
        map.get(97).add("$command_part = $cmd_param");
        map.get(88).add("$_gen11 = $cmd_param_kv $_gen12");
        map.get(90).add("$_gen12 = $cmd_param_kv $_gen12");
        map.get(90).add("$_gen12 = :_empty");
        map.get(88).add("$_gen11 = :_empty");
        map.get(59).add("$_gen13 = :string");
        map.get(59).add("$_gen13 = :_empty");
        map.get(60).add("$_gen14 = $type_e");
        map.get(60).add("$_gen14 = :_empty");
        map.get(125).add("$_gen15 = $postfix_quantifier");
        map.get(125).add("$_gen15 = :_empty");
        map.get(67).add("$cmd_param = :cmd_param_start $_gen11 $_gen13 $_gen14 :identifier $_gen15 :cmd_param_end -> CommandParameter( name=$4, type=$3, prefix=$2, attributes=$1, postfix=$5 )");
        map.get(76).add("$cmd_param_kv = :cmd_attr_hint :identifier :equal $e -> CommandParameterAttr( key=$1, value=$3 )");
        map.get(105).add("$postfix_quantifier = :qmark");
        map.get(105).add("$postfix_quantifier = :plus");
        map.get(105).add("$postfix_quantifier = :asterisk");
        map.get(120).add("$_gen16 = $output_kv $_gen17");
        map.get(84).add("$_gen17 = $output_kv $_gen17");
        map.get(84).add("$_gen17 = :_empty");
        map.get(120).add("$_gen16 = :_empty");
        map.get(109).add("$outputs = :output :lbrace $_gen16 :rbrace -> Outputs( attributes=$2 )");
        map.get(112).add("$output_kv = $type_e :identifier :equal $e -> Output( type=$0, var=$1, expression=$3 )");
        map.get(101).add("$runtime = :runtime $map -> Runtime( map=$1 )");
        map.get(58).add("$parameter_meta = :parameter_meta $map -> ParameterMeta( map=$1 )");
        map.get(85).add("$meta = :meta $map -> Meta( map=$1 )");
        map.get(70).add("$_gen18 = $kv $_gen19");
        map.get(62).add("$_gen19 = $kv $_gen19");
        map.get(62).add("$_gen19 = :_empty");
        map.get(70).add("$_gen18 = :_empty");
        map.get(79).add("$map = :lbrace $_gen18 :rbrace -> $1");
        map.get(116).add("$kv = :identifier :colon $e -> RuntimeAttribute( key=$0, value=$2 )");
        map.get(65).add("$_gen20 = $wf_body_element $_gen21");
        map.get(91).add("$_gen21 = $wf_body_element $_gen21");
        map.get(91).add("$_gen21 = :_empty");
        map.get(65).add("$_gen20 = :_empty");
        map.get(110).add("$workflow = :workflow :identifier :lbrace $_gen20 :rbrace -> Workflow( name=$1, body=$3 )");
        map.get(89).add("$wf_body_element = $call");
        map.get(89).add("$wf_body_element = $declaration");
        map.get(89).add("$wf_body_element = $while_loop");
        map.get(89).add("$wf_body_element = $if_stmt");
        map.get(89).add("$wf_body_element = $scatter");
        map.get(96).add("$_gen22 = $alias");
        map.get(96).add("$_gen22 = :_empty");
        map.get(78).add("$_gen23 = $call_body");
        map.get(78).add("$_gen23 = :_empty");
        map.get(68).add("$call = :call :ns_identifier $_gen22 $_gen23 -> Call( task=$1, alias=$2, body=$3 )");
        map.get(106).add("$_gen24 = $call_input $_gen25");
        map.get(66).add("$_gen25 = $call_input $_gen25");
        map.get(66).add("$_gen25 = :_empty");
        map.get(106).add("$_gen24 = :_empty");
        map.get(113).add("$call_body = :lbrace $_gen5 $_gen24 :rbrace -> CallBody( declarations=$1, io=$2 )");
        map.get(73).add("$_gen26 = $mapping $_gen27");
        map.get(86).add("$_gen27 = :comma $mapping $_gen27");
        map.get(86).add("$_gen27 = :_empty");
        map.get(73).add("$_gen26 = :_empty");
        map.get(72).add("$call_input = :input :colon $_gen26 -> Inputs( map=$2 )");
        map.get(114).add("$mapping = :identifier :equal $e -> IOMapping( key=$0, value=$2 )");
        map.get(115).add("$alias = :as :identifier -> $1");
        map.get(118).add("$while_loop = :while :lparen $e :rparen :lbrace $_gen20 :rbrace -> WhileLoop( expression=$2, body=$5 )");
        map.get(121).add("$if_stmt = :if :lparen $e :rparen :lbrace $_gen20 :rbrace -> If( expression=$2, body=$5 )");
        map.get(80).add("$scatter = :scatter :lparen :identifier :in $e :rparen :lbrace $_gen20 :rbrace -> Scatter( item=$2, collection=$4, body=$7 )");
        map.get(99).add("$_gen28 = $setter");
        map.get(99).add("$_gen28 = :_empty");
        map.get(92).add("$declaration = $type_e :identifier $_gen28 -> Declaration( type=$0, name=$1, expression=$2 )");
        map.get(64).add("$setter = :equal $e -> $1");
        map.get(122).add("$object_kv = :identifier :colon $e -> ObjectKV( key=$0, value=$2 )");
        map.get(102).add("$_gen29 = $type_e $_gen30");
        map.get(87).add("$_gen30 = :comma $type_e $_gen30");
        map.get(87).add("$_gen30 = :_empty");
        map.get(102).add("$_gen29 = :_empty");
        map.get(82).add("$type_e = :type <=> :lsquare $_gen29 :rsquare -> Type( name=$0, subtype=$2 )");
        map.get(82).add("$type_e = :type");
        map.get(77).add("$e = $e :double_pipe $e -> LogicalOr( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :double_ampersand $e -> LogicalAnd( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :double_equal $e -> Equals( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :not_equal $e -> NotEquals( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :lt $e -> LessThan( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :lteq $e -> LessThanOrEqual( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :gt $e -> GreaterThan( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :gteq $e -> GreaterThanOrEqual( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :plus $e -> Add( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :dash $e -> Subtract( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :asterisk $e -> Multiply( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :slash $e -> Divide( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = $e :percent $e -> Remainder( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = :not $e -> LogicalNot( expression=$1 )");
        map.get(77).add("$e = :plus $e -> UnaryPlus( expression=$1 )");
        map.get(77).add("$e = :dash $e -> UnaryNegation( expression=$1 )");
        map.get(123).add("$_gen31 = $e $_gen32");
        map.get(107).add("$_gen32 = :comma $e $_gen32");
        map.get(107).add("$_gen32 = :_empty");
        map.get(123).add("$_gen31 = :_empty");
        map.get(77).add("$e = :identifier <=> :lparen $_gen31 :rparen -> FunctionCall( name=$0, params=$2 )");
        map.get(77).add("$e = :identifier <=> :lsquare $e :rsquare -> ArrayIndex( lhs=$0, rhs=$2 )");
        map.get(77).add("$e = :identifier <=> :dot :identifier -> MemberAccess( lhs=$0, rhs=$2 )");
        map.get(103).add("$_gen33 = $object_kv $_gen34");
        map.get(75).add("$_gen34 = :comma $object_kv $_gen34");
        map.get(75).add("$_gen34 = :_empty");
        map.get(103).add("$_gen33 = :_empty");
        map.get(77).add("$e = :object :lbrace $_gen33 :rbrace -> ObjectLiteral( map=$2 )");
        map.get(77).add("$e = :lsquare $_gen31 :rsquare -> ArrayLiteral( values=$1 )");
        map.get(77).add("$e = :lparen $e :rparen -> $1");
        map.get(77).add("$e = :string");
        map.get(77).add("$e = :identifier");
        map.get(77).add("$e = :boolean");
        map.get(77).add("$e = :integer");
        map.get(77).add("$e = :float");
        map.get(77).add("$e = :dquote_string");
        map.get(77).add("$e = :squote_string");
        nonterminal_rules = Collections.unmodifiableMap(map);
    }
    static {
        Map<Integer, String> map = new HashMap<Integer, String>();
        map.put(new Integer(0), "$_gen0 = $import $_gen1");
        map.put(new Integer(1), "$_gen1 = $import $_gen1");
        map.put(new Integer(2), "$_gen1 = :_empty");
        map.put(new Integer(3), "$_gen0 = :_empty");
        map.put(new Integer(4), "$_gen2 = $workflow_or_task $_gen3");
        map.put(new Integer(5), "$_gen3 = $workflow_or_task $_gen3");
        map.put(new Integer(6), "$_gen3 = :_empty");
        map.put(new Integer(7), "$_gen2 = :_empty");
        map.put(new Integer(8), "$document = $_gen0 $_gen2 -> Document( imports=$0, definitions=$1 )");
        map.put(new Integer(9), "$workflow_or_task = $workflow");
        map.put(new Integer(10), "$workflow_or_task = $task");
        map.put(new Integer(11), "$_gen4 = $import_namespace");
        map.put(new Integer(12), "$_gen4 = :_empty");
        map.put(new Integer(13), "$import = :import :string $_gen4 -> Import( uri=$1, namespace=$2 )");
        map.put(new Integer(14), "$import_namespace = :as :identifier -> $1");
        map.put(new Integer(15), "$_gen5 = $declaration $_gen6");
        map.put(new Integer(16), "$_gen6 = $declaration $_gen6");
        map.put(new Integer(17), "$_gen6 = :_empty");
        map.put(new Integer(18), "$_gen5 = :_empty");
        map.put(new Integer(19), "$_gen7 = $sections $_gen8");
        map.put(new Integer(20), "$_gen8 = $sections $_gen8");
        map.put(new Integer(21), "$_gen8 = :_empty");
        map.put(new Integer(22), "$_gen7 = :_empty");
        map.put(new Integer(23), "$task = :task :identifier :lbrace $_gen5 $_gen7 :rbrace -> Task( name=$1, declarations=$3, sections=$4 )");
        map.put(new Integer(24), "$sections = $command");
        map.put(new Integer(25), "$sections = $outputs");
        map.put(new Integer(26), "$sections = $runtime");
        map.put(new Integer(27), "$sections = $parameter_meta");
        map.put(new Integer(28), "$sections = $meta");
        map.put(new Integer(29), "$_gen9 = $command_part $_gen10");
        map.put(new Integer(30), "$_gen10 = $command_part $_gen10");
        map.put(new Integer(31), "$_gen10 = :_empty");
        map.put(new Integer(32), "$_gen9 = :_empty");
        map.put(new Integer(33), "$command = :raw_command :raw_cmd_start $_gen9 :raw_cmd_end -> RawCommand( parts=$2 )");
        map.put(new Integer(34), "$command_part = :cmd_part");
        map.put(new Integer(35), "$command_part = $cmd_param");
        map.put(new Integer(36), "$_gen11 = $cmd_param_kv $_gen12");
        map.put(new Integer(37), "$_gen12 = $cmd_param_kv $_gen12");
        map.put(new Integer(38), "$_gen12 = :_empty");
        map.put(new Integer(39), "$_gen11 = :_empty");
        map.put(new Integer(40), "$_gen13 = :string");
        map.put(new Integer(41), "$_gen13 = :_empty");
        map.put(new Integer(42), "$_gen14 = $type_e");
        map.put(new Integer(43), "$_gen14 = :_empty");
        map.put(new Integer(44), "$_gen15 = $postfix_quantifier");
        map.put(new Integer(45), "$_gen15 = :_empty");
        map.put(new Integer(46), "$cmd_param = :cmd_param_start $_gen11 $_gen13 $_gen14 :identifier $_gen15 :cmd_param_end -> CommandParameter( name=$4, type=$3, prefix=$2, attributes=$1, postfix=$5 )");
        map.put(new Integer(47), "$cmd_param_kv = :cmd_attr_hint :identifier :equal $e -> CommandParameterAttr( key=$1, value=$3 )");
        map.put(new Integer(48), "$postfix_quantifier = :qmark");
        map.put(new Integer(49), "$postfix_quantifier = :plus");
        map.put(new Integer(50), "$postfix_quantifier = :asterisk");
        map.put(new Integer(51), "$_gen16 = $output_kv $_gen17");
        map.put(new Integer(52), "$_gen17 = $output_kv $_gen17");
        map.put(new Integer(53), "$_gen17 = :_empty");
        map.put(new Integer(54), "$_gen16 = :_empty");
        map.put(new Integer(55), "$outputs = :output :lbrace $_gen16 :rbrace -> Outputs( attributes=$2 )");
        map.put(new Integer(56), "$output_kv = $type_e :identifier :equal $e -> Output( type=$0, var=$1, expression=$3 )");
        map.put(new Integer(57), "$runtime = :runtime $map -> Runtime( map=$1 )");
        map.put(new Integer(58), "$parameter_meta = :parameter_meta $map -> ParameterMeta( map=$1 )");
        map.put(new Integer(59), "$meta = :meta $map -> Meta( map=$1 )");
        map.put(new Integer(60), "$_gen18 = $kv $_gen19");
        map.put(new Integer(61), "$_gen19 = $kv $_gen19");
        map.put(new Integer(62), "$_gen19 = :_empty");
        map.put(new Integer(63), "$_gen18 = :_empty");
        map.put(new Integer(64), "$map = :lbrace $_gen18 :rbrace -> $1");
        map.put(new Integer(65), "$kv = :identifier :colon $e -> RuntimeAttribute( key=$0, value=$2 )");
        map.put(new Integer(66), "$_gen20 = $wf_body_element $_gen21");
        map.put(new Integer(67), "$_gen21 = $wf_body_element $_gen21");
        map.put(new Integer(68), "$_gen21 = :_empty");
        map.put(new Integer(69), "$_gen20 = :_empty");
        map.put(new Integer(70), "$workflow = :workflow :identifier :lbrace $_gen20 :rbrace -> Workflow( name=$1, body=$3 )");
        map.put(new Integer(71), "$wf_body_element = $call");
        map.put(new Integer(72), "$wf_body_element = $declaration");
        map.put(new Integer(73), "$wf_body_element = $while_loop");
        map.put(new Integer(74), "$wf_body_element = $if_stmt");
        map.put(new Integer(75), "$wf_body_element = $scatter");
        map.put(new Integer(76), "$_gen22 = $alias");
        map.put(new Integer(77), "$_gen22 = :_empty");
        map.put(new Integer(78), "$_gen23 = $call_body");
        map.put(new Integer(79), "$_gen23 = :_empty");
        map.put(new Integer(80), "$call = :call :ns_identifier $_gen22 $_gen23 -> Call( task=$1, alias=$2, body=$3 )");
        map.put(new Integer(81), "$_gen24 = $call_input $_gen25");
        map.put(new Integer(82), "$_gen25 = $call_input $_gen25");
        map.put(new Integer(83), "$_gen25 = :_empty");
        map.put(new Integer(84), "$_gen24 = :_empty");
        map.put(new Integer(85), "$call_body = :lbrace $_gen5 $_gen24 :rbrace -> CallBody( declarations=$1, io=$2 )");
        map.put(new Integer(86), "$_gen26 = $mapping $_gen27");
        map.put(new Integer(87), "$_gen27 = :comma $mapping $_gen27");
        map.put(new Integer(88), "$_gen27 = :_empty");
        map.put(new Integer(89), "$_gen26 = :_empty");
        map.put(new Integer(90), "$call_input = :input :colon $_gen26 -> Inputs( map=$2 )");
        map.put(new Integer(91), "$mapping = :identifier :equal $e -> IOMapping( key=$0, value=$2 )");
        map.put(new Integer(92), "$alias = :as :identifier -> $1");
        map.put(new Integer(93), "$while_loop = :while :lparen $e :rparen :lbrace $_gen20 :rbrace -> WhileLoop( expression=$2, body=$5 )");
        map.put(new Integer(94), "$if_stmt = :if :lparen $e :rparen :lbrace $_gen20 :rbrace -> If( expression=$2, body=$5 )");
        map.put(new Integer(95), "$scatter = :scatter :lparen :identifier :in $e :rparen :lbrace $_gen20 :rbrace -> Scatter( item=$2, collection=$4, body=$7 )");
        map.put(new Integer(96), "$_gen28 = $setter");
        map.put(new Integer(97), "$_gen28 = :_empty");
        map.put(new Integer(98), "$declaration = $type_e :identifier $_gen28 -> Declaration( type=$0, name=$1, expression=$2 )");
        map.put(new Integer(99), "$setter = :equal $e -> $1");
        map.put(new Integer(100), "$object_kv = :identifier :colon $e -> ObjectKV( key=$0, value=$2 )");
        map.put(new Integer(101), "$_gen29 = $type_e $_gen30");
        map.put(new Integer(102), "$_gen30 = :comma $type_e $_gen30");
        map.put(new Integer(103), "$_gen30 = :_empty");
        map.put(new Integer(104), "$_gen29 = :_empty");
        map.put(new Integer(105), "$type_e = :type <=> :lsquare $_gen29 :rsquare -> Type( name=$0, subtype=$2 )");
        map.put(new Integer(106), "$type_e = :type");
        map.put(new Integer(107), "$e = $e :double_pipe $e -> LogicalOr( lhs=$0, rhs=$2 )");
        map.put(new Integer(108), "$e = $e :double_ampersand $e -> LogicalAnd( lhs=$0, rhs=$2 )");
        map.put(new Integer(109), "$e = $e :double_equal $e -> Equals( lhs=$0, rhs=$2 )");
        map.put(new Integer(110), "$e = $e :not_equal $e -> NotEquals( lhs=$0, rhs=$2 )");
        map.put(new Integer(111), "$e = $e :lt $e -> LessThan( lhs=$0, rhs=$2 )");
        map.put(new Integer(112), "$e = $e :lteq $e -> LessThanOrEqual( lhs=$0, rhs=$2 )");
        map.put(new Integer(113), "$e = $e :gt $e -> GreaterThan( lhs=$0, rhs=$2 )");
        map.put(new Integer(114), "$e = $e :gteq $e -> GreaterThanOrEqual( lhs=$0, rhs=$2 )");
        map.put(new Integer(115), "$e = $e :plus $e -> Add( lhs=$0, rhs=$2 )");
        map.put(new Integer(116), "$e = $e :dash $e -> Subtract( lhs=$0, rhs=$2 )");
        map.put(new Integer(117), "$e = $e :asterisk $e -> Multiply( lhs=$0, rhs=$2 )");
        map.put(new Integer(118), "$e = $e :slash $e -> Divide( lhs=$0, rhs=$2 )");
        map.put(new Integer(119), "$e = $e :percent $e -> Remainder( lhs=$0, rhs=$2 )");
        map.put(new Integer(120), "$e = :not $e -> LogicalNot( expression=$1 )");
        map.put(new Integer(121), "$e = :plus $e -> UnaryPlus( expression=$1 )");
        map.put(new Integer(122), "$e = :dash $e -> UnaryNegation( expression=$1 )");
        map.put(new Integer(123), "$_gen31 = $e $_gen32");
        map.put(new Integer(124), "$_gen32 = :comma $e $_gen32");
        map.put(new Integer(125), "$_gen32 = :_empty");
        map.put(new Integer(126), "$_gen31 = :_empty");
        map.put(new Integer(127), "$e = :identifier <=> :lparen $_gen31 :rparen -> FunctionCall( name=$0, params=$2 )");
        map.put(new Integer(128), "$e = :identifier <=> :lsquare $e :rsquare -> ArrayIndex( lhs=$0, rhs=$2 )");
        map.put(new Integer(129), "$e = :identifier <=> :dot :identifier -> MemberAccess( lhs=$0, rhs=$2 )");
        map.put(new Integer(130), "$_gen33 = $object_kv $_gen34");
        map.put(new Integer(131), "$_gen34 = :comma $object_kv $_gen34");
        map.put(new Integer(132), "$_gen34 = :_empty");
        map.put(new Integer(133), "$_gen33 = :_empty");
        map.put(new Integer(134), "$e = :object :lbrace $_gen33 :rbrace -> ObjectLiteral( map=$2 )");
        map.put(new Integer(135), "$e = :lsquare $_gen31 :rsquare -> ArrayLiteral( values=$1 )");
        map.put(new Integer(136), "$e = :lparen $e :rparen -> $1");
        map.put(new Integer(137), "$e = :string");
        map.put(new Integer(138), "$e = :identifier");
        map.put(new Integer(139), "$e = :boolean");
        map.put(new Integer(140), "$e = :integer");
        map.put(new Integer(141), "$e = :float");
        map.put(new Integer(142), "$e = :dquote_string");
        map.put(new Integer(143), "$e = :squote_string");
        rules = Collections.unmodifiableMap(map);
    }
    public static boolean is_terminal(int id) {
        return 0 <= id && id <= 57;
    }
    public ParseTree parse(TokenStream tokens) throws SyntaxError {
        return parse(tokens, new DefaultSyntaxErrorFormatter());
    }
    public ParseTree parse(List<Terminal> tokens) throws SyntaxError {
        return parse(new TokenStream(tokens));
    }
    public ParseTree parse(TokenStream tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(tokens, error_formatter);
        ParseTree tree = parse_document(ctx);
        if (ctx.tokens.current() != null) {
            StackTraceElement[] stack = Thread.currentThread().getStackTrace();
            throw new SyntaxError(ctx.error_formatter.excessTokens(stack[1].getMethodName(), ctx.tokens.current()));
        }
        return tree;
    }
    public ParseTree parse(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        return parse(new TokenStream(tokens), error_formatter);
    }
    private static Terminal expect(ParserContext ctx, TerminalIdentifier expecting) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.noMoreTokens(ctx.nonterminal, expecting, ctx.tokens.last()));
        }
        if (current.getId() != expecting.id()) {
            ArrayList<TerminalIdentifier> expectedList = new ArrayList<TerminalIdentifier>();
            expectedList.add(expecting);
            throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(ctx.nonterminal, current, expectedList, ctx.rule));
        }
        Terminal next = ctx.tokens.advance();
        if ( next != null && !is_terminal(next.getId()) ) {
            throw new SyntaxError(ctx.error_formatter.invalidTerminal(ctx.nonterminal, next));
        }
        return current;
    }
    private static Map<Integer, Integer> infix_binding_power_e;
    private static Map<Integer, Integer> prefix_binding_power_e;
    static {
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        map.put(44, 2000); /* $e = $e :double_pipe $e -> LogicalOr( lhs=$0, rhs=$2 ) */
        map.put(43, 3000); /* $e = $e :double_ampersand $e -> LogicalAnd( lhs=$0, rhs=$2 ) */
        map.put(23, 4000); /* $e = $e :double_equal $e -> Equals( lhs=$0, rhs=$2 ) */
        map.put(33, 4000); /* $e = $e :not_equal $e -> NotEquals( lhs=$0, rhs=$2 ) */
        map.put(15, 5000); /* $e = $e :lt $e -> LessThan( lhs=$0, rhs=$2 ) */
        map.put(57, 5000); /* $e = $e :lteq $e -> LessThanOrEqual( lhs=$0, rhs=$2 ) */
        map.put(9, 5000); /* $e = $e :gt $e -> GreaterThan( lhs=$0, rhs=$2 ) */
        map.put(54, 5000); /* $e = $e :gteq $e -> GreaterThanOrEqual( lhs=$0, rhs=$2 ) */
        map.put(29, 6000); /* $e = $e :plus $e -> Add( lhs=$0, rhs=$2 ) */
        map.put(38, 6000); /* $e = $e :dash $e -> Subtract( lhs=$0, rhs=$2 ) */
        map.put(8, 7000); /* $e = $e :asterisk $e -> Multiply( lhs=$0, rhs=$2 ) */
        map.put(22, 7000); /* $e = $e :slash $e -> Divide( lhs=$0, rhs=$2 ) */
        map.put(56, 7000); /* $e = $e :percent $e -> Remainder( lhs=$0, rhs=$2 ) */
        map.put(20, 9000); /* $e = :identifier <=> :lparen list(nt=$e, sep=:comma, min=0, sep_terminates=False) :rparen -> FunctionCall( name=$0, params=$2 ) */
        map.put(28, 10000); /* $e = :identifier <=> :lsquare $e :rsquare -> ArrayIndex( lhs=$0, rhs=$2 ) */
        map.put(5, 11000); /* $e = :identifier <=> :dot :identifier -> MemberAccess( lhs=$0, rhs=$2 ) */
        infix_binding_power_e = Collections.unmodifiableMap(map);
    }
    static {
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        map.put(16, 8000); /* $e = :not $e -> LogicalNot( expression=$1 ) */
        map.put(29, 8000); /* $e = :plus $e -> UnaryPlus( expression=$1 ) */
        map.put(38, 8000); /* $e = :dash $e -> UnaryNegation( expression=$1 ) */
        prefix_binding_power_e = Collections.unmodifiableMap(map);
    }
    static int get_infix_binding_power_e(int terminal_id) {
        if (infix_binding_power_e.containsKey(terminal_id)) {
            return infix_binding_power_e.get(terminal_id);
        }
        return 0;
    }
    static int get_prefix_binding_power_e(int terminal_id) {
        if (prefix_binding_power_e.containsKey(terminal_id)) {
            return prefix_binding_power_e.get(terminal_id);
        }
        return 0;
    }
    public ParseTree parse_e(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_e_internal(ctx, 0);
    }
    public static ParseTree parse_e(ParserContext ctx) throws SyntaxError {
        return parse_e_internal(ctx, 0);
    }
    public static ParseTree parse_e_internal(ParserContext ctx, int rbp) throws SyntaxError {
        ParseTree left = nud_e(ctx);
        if ( left instanceof ParseTree ) {
            left.setExpr(true);
            left.setNud(true);
        }
        while (ctx.tokens.current() != null && rbp < get_infix_binding_power_e(ctx.tokens.current().getId())) {
            left = led_e(left, ctx);
        }
        if (left != null) {
            left.setExpr(true);
        }
        return left;
    }
    private static ParseTree nud_e(ParserContext ctx) throws SyntaxError {
        ParseTree tree = new ParseTree( new NonTerminal(77, "e") );
        Terminal current = ctx.tokens.current();
        ctx.nonterminal = "e";
        if (current == null) {
            return tree;
        }
        else if (rule_first.get(120).contains(terminal_map.get(current.getId()))) {
            /* (120) $e = :not $e -> LogicalNot( expression=$1 ) */
            ctx.rule = rules.get(120);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("expression", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("LogicalNot", parameters));
            tree.setNudMorphemeCount(2);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_NOT));
            tree.add(parse_e_internal(ctx, get_prefix_binding_power_e(16)));
            tree.setPrefix(true);
        }
        else if (rule_first.get(121).contains(terminal_map.get(current.getId()))) {
            /* (121) $e = :plus $e -> UnaryPlus( expression=$1 ) */
            ctx.rule = rules.get(121);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("expression", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("UnaryPlus", parameters));
            tree.setNudMorphemeCount(2);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_PLUS));
            tree.add(parse_e_internal(ctx, get_prefix_binding_power_e(29)));
            tree.setPrefix(true);
        }
        else if (rule_first.get(122).contains(terminal_map.get(current.getId()))) {
            /* (122) $e = :dash $e -> UnaryNegation( expression=$1 ) */
            ctx.rule = rules.get(122);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("expression", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("UnaryNegation", parameters));
            tree.setNudMorphemeCount(2);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DASH));
            tree.add(parse_e_internal(ctx, get_prefix_binding_power_e(38)));
            tree.setPrefix(true);
        }
        else if (rule_first.get(127).contains(terminal_map.get(current.getId()))) {
            /* (127) $e = :identifier <=> :lparen $_gen31 :rparen -> FunctionCall( name=$0, params=$2 ) */
            ctx.rule = rules.get(127);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER));
        }
        else if (rule_first.get(128).contains(terminal_map.get(current.getId()))) {
            /* (128) $e = :identifier <=> :lsquare $e :rsquare -> ArrayIndex( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(128);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER));
        }
        else if (rule_first.get(129).contains(terminal_map.get(current.getId()))) {
            /* (129) $e = :identifier <=> :dot :identifier -> MemberAccess( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(129);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER));
        }
        else if (rule_first.get(134).contains(terminal_map.get(current.getId()))) {
            /* (134) $e = :object :lbrace $_gen33 :rbrace -> ObjectLiteral( map=$2 ) */
            ctx.rule = rules.get(134);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("map", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("ObjectLiteral", parameters));
            tree.setNudMorphemeCount(4);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_OBJECT));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE));
            tree.add(parse__gen33(ctx));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE));
        }
        else if (rule_first.get(135).contains(terminal_map.get(current.getId()))) {
            /* (135) $e = :lsquare $_gen31 :rsquare -> ArrayLiteral( values=$1 ) */
            ctx.rule = rules.get(135);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("values", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("ArrayLiteral", parameters));
            tree.setNudMorphemeCount(3);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LSQUARE));
            tree.add(parse__gen31(ctx));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_RSQUARE));
        }
        else if (rule_first.get(136).contains(terminal_map.get(current.getId()))) {
            /* (136) $e = :lparen $e :rparen -> $1 */
            ctx.rule = rules.get(136);
            tree.setAstTransformation(new AstTransformSubstitution(1));
            tree.setNudMorphemeCount(3);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LPAREN));
            tree.add(parse_e(ctx));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_RPAREN));
        }
        else if (rule_first.get(137).contains(terminal_map.get(current.getId()))) {
            /* (137) $e = :string */
            ctx.rule = rules.get(137);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_STRING));
        }
        else if (rule_first.get(138).contains(terminal_map.get(current.getId()))) {
            /* (138) $e = :identifier */
            ctx.rule = rules.get(138);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER));
        }
        else if (rule_first.get(139).contains(terminal_map.get(current.getId()))) {
            /* (139) $e = :boolean */
            ctx.rule = rules.get(139);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_BOOLEAN));
        }
        else if (rule_first.get(140).contains(terminal_map.get(current.getId()))) {
            /* (140) $e = :integer */
            ctx.rule = rules.get(140);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_INTEGER));
        }
        else if (rule_first.get(141).contains(terminal_map.get(current.getId()))) {
            /* (141) $e = :float */
            ctx.rule = rules.get(141);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_FLOAT));
        }
        else if (rule_first.get(142).contains(terminal_map.get(current.getId()))) {
            /* (142) $e = :dquote_string */
            ctx.rule = rules.get(142);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DQUOTE_STRING));
        }
        else if (rule_first.get(143).contains(terminal_map.get(current.getId()))) {
            /* (143) $e = :squote_string */
            ctx.rule = rules.get(143);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_SQUOTE_STRING));
        }
        return tree;
    }
    private static ParseTree led_e(ParseTree left, ParserContext ctx) throws SyntaxError {
        ParseTree tree = new ParseTree( new NonTerminal(77, "e") );
        Terminal current = ctx.tokens.current();
        ctx.nonterminal = "e";
        int modifier;
        if (current.getId() == 44) {
            /* $e = $e :double_pipe $e -> LogicalOr( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(107);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("LogicalOr", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DOUBLE_PIPE));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(44) - modifier));
            return tree;
        }
        if (current.getId() == 43) {
            /* $e = $e :double_ampersand $e -> LogicalAnd( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(108);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("LogicalAnd", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DOUBLE_AMPERSAND));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(43) - modifier));
            return tree;
        }
        if (current.getId() == 23) {
            /* $e = $e :double_equal $e -> Equals( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(109);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Equals", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DOUBLE_EQUAL));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(23) - modifier));
            return tree;
        }
        if (current.getId() == 33) {
            /* $e = $e :not_equal $e -> NotEquals( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(110);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("NotEquals", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_NOT_EQUAL));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(33) - modifier));
            return tree;
        }
        if (current.getId() == 15) {
            /* $e = $e :lt $e -> LessThan( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(111);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("LessThan", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LT));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(15) - modifier));
            return tree;
        }
        if (current.getId() == 57) {
            /* $e = $e :lteq $e -> LessThanOrEqual( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(112);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("LessThanOrEqual", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LTEQ));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(57) - modifier));
            return tree;
        }
        if (current.getId() == 9) {
            /* $e = $e :gt $e -> GreaterThan( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(113);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("GreaterThan", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_GT));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(9) - modifier));
            return tree;
        }
        if (current.getId() == 54) {
            /* $e = $e :gteq $e -> GreaterThanOrEqual( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(114);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("GreaterThanOrEqual", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_GTEQ));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(54) - modifier));
            return tree;
        }
        if (current.getId() == 29) {
            /* $e = $e :plus $e -> Add( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(115);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Add", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_PLUS));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(29) - modifier));
            return tree;
        }
        if (current.getId() == 38) {
            /* $e = $e :dash $e -> Subtract( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(116);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Subtract", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DASH));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(38) - modifier));
            return tree;
        }
        if (current.getId() == 8) {
            /* $e = $e :asterisk $e -> Multiply( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(117);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Multiply", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_ASTERISK));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(8) - modifier));
            return tree;
        }
        if (current.getId() == 22) {
            /* $e = $e :slash $e -> Divide( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(118);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Divide", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_SLASH));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(22) - modifier));
            return tree;
        }
        if (current.getId() == 56) {
            /* $e = $e :percent $e -> Remainder( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(119);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Remainder", parameters));
            tree.setExprNud(true);
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_PERCENT));
            modifier = 0;
            tree.setInfix(true);
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(56) - modifier));
            return tree;
        }
        if (current.getId() == 20) {
            /* $e = :identifier <=> :lparen $_gen31 :rparen -> FunctionCall( name=$0, params=$2 ) */
            ctx.rule = rules.get(127);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("name", 0);
            parameters.put("params", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("FunctionCall", parameters));
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LPAREN));
            tree.add(parse__gen31(ctx));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_RPAREN));
            return tree;
        }
        if (current.getId() == 28) {
            /* $e = :identifier <=> :lsquare $e :rsquare -> ArrayIndex( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(128);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("ArrayIndex", parameters));
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LSQUARE));
            modifier = 0;
            tree.add(parse_e_internal(ctx, get_infix_binding_power_e(28) - modifier));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_RSQUARE));
            return tree;
        }
        if (current.getId() == 5) {
            /* $e = :identifier <=> :dot :identifier -> MemberAccess( lhs=$0, rhs=$2 ) */
            ctx.rule = rules.get(129);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("lhs", 0);
            parameters.put("rhs", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("MemberAccess", parameters));
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_DOT));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER));
            return tree;
        }
        return tree;
    }
    private static Map<Integer, Integer> infix_binding_power_type_e;
    private static Map<Integer, Integer> prefix_binding_power_type_e;
    static {
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        map.put(28, 1000); /* $type_e = :type <=> :lsquare list(nt=$type_e, sep=:comma, min=0, sep_terminates=False) :rsquare -> Type( name=$0, subtype=$2 ) */
        infix_binding_power_type_e = Collections.unmodifiableMap(map);
    }
    static {
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        prefix_binding_power_type_e = Collections.unmodifiableMap(map);
    }
    static int get_infix_binding_power_type_e(int terminal_id) {
        if (infix_binding_power_type_e.containsKey(terminal_id)) {
            return infix_binding_power_type_e.get(terminal_id);
        }
        return 0;
    }
    static int get_prefix_binding_power_type_e(int terminal_id) {
        if (prefix_binding_power_type_e.containsKey(terminal_id)) {
            return prefix_binding_power_type_e.get(terminal_id);
        }
        return 0;
    }
    public ParseTree parse_type_e(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_type_e_internal(ctx, 0);
    }
    public static ParseTree parse_type_e(ParserContext ctx) throws SyntaxError {
        return parse_type_e_internal(ctx, 0);
    }
    public static ParseTree parse_type_e_internal(ParserContext ctx, int rbp) throws SyntaxError {
        ParseTree left = nud_type_e(ctx);
        if ( left instanceof ParseTree ) {
            left.setExpr(true);
            left.setNud(true);
        }
        while (ctx.tokens.current() != null && rbp < get_infix_binding_power_type_e(ctx.tokens.current().getId())) {
            left = led_type_e(left, ctx);
        }
        if (left != null) {
            left.setExpr(true);
        }
        return left;
    }
    private static ParseTree nud_type_e(ParserContext ctx) throws SyntaxError {
        ParseTree tree = new ParseTree( new NonTerminal(82, "type_e") );
        Terminal current = ctx.tokens.current();
        ctx.nonterminal = "type_e";
        if (current == null) {
            return tree;
        }
        if (rule_first.get(105).contains(terminal_map.get(current.getId()))) {
            /* (105) $type_e = :type <=> :lsquare $_gen29 :rsquare -> Type( name=$0, subtype=$2 ) */
            ctx.rule = rules.get(105);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_TYPE));
        }
        else if (rule_first.get(106).contains(terminal_map.get(current.getId()))) {
            /* (106) $type_e = :type */
            ctx.rule = rules.get(106);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            tree.setNudMorphemeCount(1);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_TYPE));
        }
        return tree;
    }
    private static ParseTree led_type_e(ParseTree left, ParserContext ctx) throws SyntaxError {
        ParseTree tree = new ParseTree( new NonTerminal(82, "type_e") );
        Terminal current = ctx.tokens.current();
        ctx.nonterminal = "type_e";
        int modifier;
        if (current.getId() == 28) {
            /* $type_e = :type <=> :lsquare $_gen29 :rsquare -> Type( name=$0, subtype=$2 ) */
            ctx.rule = rules.get(105);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("name", 0);
            parameters.put("subtype", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Type", parameters));
            tree.add(left);
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_LSQUARE));
            tree.add(parse__gen29(ctx));
            tree.add(expect(ctx, WdlTerminalIdentifier.TERMINAL_RSQUARE));
            return tree;
        }
        return tree;
    }
    public ParseTree parse_parameter_meta(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_parameter_meta(ctx);
    }
    private static ParseTree parse_parameter_meta(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[0][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(58, "parameter_meta"));
        ctx.nonterminal = "parameter_meta";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "parameter_meta",
                nonterminal_first.get(58),
                nonterminal_rules.get(58)
            ));
        }
        if (rule == 58) {
            /* $parameter_meta = :parameter_meta $map -> ParameterMeta( map=$1 ) */
            ctx.rule = rules.get(58);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("map", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("ParameterMeta", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_PARAMETER_META);
            tree.add(next);
            subtree = parse_map(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "parameter_meta",
            current,
            nonterminal_first.get(58),
            rules.get(58)
        ));
    }
    public ParseTree parse__gen13(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen13(ctx);
    }
    private static ParseTree parse__gen13(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[1][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(59, "_gen13"));
        ctx.nonterminal = "_gen13";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(59).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(59).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 40) {
            /* $_gen13 = :string */
            ctx.rule = rules.get(40);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_STRING);
            tree.add(next);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen14(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen14(ctx);
    }
    private static ParseTree parse__gen14(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[2][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(60, "_gen14"));
        ctx.nonterminal = "_gen14";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(60).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(60).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 42) {
            /* $_gen14 = $type_e */
            ctx.rule = rules.get(42);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_type_e(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen7(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen7(ctx);
    }
    private static ParseTree parse__gen7(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[3][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(61, "_gen7"));
        ctx.nonterminal = "_gen7";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(61).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(61).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 19) {
            /* $_gen7 = $sections $_gen8 */
            ctx.rule = rules.get(19);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_sections(ctx);
            tree.add(subtree);
            subtree = parse__gen8(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen19(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen19(ctx);
    }
    private static ParseTree parse__gen19(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[4][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(62, "_gen19"));
        ctx.nonterminal = "_gen19";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(62).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(62).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 61) {
            /* $_gen19 = $kv $_gen19 */
            ctx.rule = rules.get(61);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen19(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen5(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen5(ctx);
    }
    private static ParseTree parse__gen5(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[5][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(63, "_gen5"));
        ctx.nonterminal = "_gen5";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(63).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(63).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 15) {
            /* $_gen5 = $declaration $_gen6 */
            ctx.rule = rules.get(15);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_declaration(ctx);
            tree.add(subtree);
            subtree = parse__gen6(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_setter(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_setter(ctx);
    }
    private static ParseTree parse_setter(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[6][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(64, "setter"));
        ctx.nonterminal = "setter";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "setter",
                nonterminal_first.get(64),
                nonterminal_rules.get(64)
            ));
        }
        if (rule == 99) {
            /* $setter = :equal $e -> $1 */
            ctx.rule = rules.get(99);
            tree.setAstTransformation(new AstTransformSubstitution(1));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_EQUAL);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "setter",
            current,
            nonterminal_first.get(64),
            rules.get(99)
        ));
    }
    public ParseTree parse__gen20(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen20(ctx);
    }
    private static ParseTree parse__gen20(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[7][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(65, "_gen20"));
        ctx.nonterminal = "_gen20";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(65).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(65).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 66) {
            /* $_gen20 = $wf_body_element $_gen21 */
            ctx.rule = rules.get(66);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_wf_body_element(ctx);
            tree.add(subtree);
            subtree = parse__gen21(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen25(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen25(ctx);
    }
    private static ParseTree parse__gen25(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[8][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(66, "_gen25"));
        ctx.nonterminal = "_gen25";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(66).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(66).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 82) {
            /* $_gen25 = $call_input $_gen25 */
            ctx.rule = rules.get(82);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_call_input(ctx);
            tree.add(subtree);
            subtree = parse__gen25(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_cmd_param(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_cmd_param(ctx);
    }
    private static ParseTree parse_cmd_param(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[9][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(67, "cmd_param"));
        ctx.nonterminal = "cmd_param";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "cmd_param",
                nonterminal_first.get(67),
                nonterminal_rules.get(67)
            ));
        }
        if (rule == 46) {
            /* $cmd_param = :cmd_param_start $_gen11 $_gen13 $_gen14 :identifier $_gen15 :cmd_param_end -> CommandParameter( name=$4, type=$3, prefix=$2, attributes=$1, postfix=$5 ) */
            ctx.rule = rules.get(46);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("name", 4);
            parameters.put("type", 3);
            parameters.put("prefix", 2);
            parameters.put("attributes", 1);
            parameters.put("postfix", 5);
            tree.setAstTransformation(new AstTransformNodeCreator("CommandParameter", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START);
            tree.add(next);
            subtree = parse__gen11(ctx);
            tree.add(subtree);
            subtree = parse__gen13(ctx);
            tree.add(subtree);
            subtree = parse__gen14(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            subtree = parse__gen15(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_CMD_PARAM_END);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "cmd_param",
            current,
            nonterminal_first.get(67),
            rules.get(46)
        ));
    }
    public ParseTree parse_call(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_call(ctx);
    }
    private static ParseTree parse_call(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[10][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(68, "call"));
        ctx.nonterminal = "call";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "call",
                nonterminal_first.get(68),
                nonterminal_rules.get(68)
            ));
        }
        if (rule == 80) {
            /* $call = :call :ns_identifier $_gen22 $_gen23 -> Call( task=$1, alias=$2, body=$3 ) */
            ctx.rule = rules.get(80);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("task", 1);
            parameters.put("alias", 2);
            parameters.put("body", 3);
            tree.setAstTransformation(new AstTransformNodeCreator("Call", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_CALL);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_NS_IDENTIFIER);
            tree.add(next);
            subtree = parse__gen22(ctx);
            tree.add(subtree);
            subtree = parse__gen23(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "call",
            current,
            nonterminal_first.get(68),
            rules.get(80)
        ));
    }
    public ParseTree parse_import(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_import(ctx);
    }
    private static ParseTree parse_import(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[11][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(69, "import"));
        ctx.nonterminal = "import";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "import",
                nonterminal_first.get(69),
                nonterminal_rules.get(69)
            ));
        }
        if (rule == 13) {
            /* $import = :import :string $_gen4 -> Import( uri=$1, namespace=$2 ) */
            ctx.rule = rules.get(13);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("uri", 1);
            parameters.put("namespace", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Import", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IMPORT);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_STRING);
            tree.add(next);
            subtree = parse__gen4(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "import",
            current,
            nonterminal_first.get(69),
            rules.get(13)
        ));
    }
    public ParseTree parse__gen18(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen18(ctx);
    }
    private static ParseTree parse__gen18(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[12][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(70, "_gen18"));
        ctx.nonterminal = "_gen18";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(70).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(70).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 60) {
            /* $_gen18 = $kv $_gen19 */
            ctx.rule = rules.get(60);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen19(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen1(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen1(ctx);
    }
    private static ParseTree parse__gen1(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[13][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(71, "_gen1"));
        ctx.nonterminal = "_gen1";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(71).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(71).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 1) {
            /* $_gen1 = $import $_gen1 */
            ctx.rule = rules.get(1);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_import(ctx);
            tree.add(subtree);
            subtree = parse__gen1(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_call_input(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_call_input(ctx);
    }
    private static ParseTree parse_call_input(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[14][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(72, "call_input"));
        ctx.nonterminal = "call_input";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "call_input",
                nonterminal_first.get(72),
                nonterminal_rules.get(72)
            ));
        }
        if (rule == 90) {
            /* $call_input = :input :colon $_gen26 -> Inputs( map=$2 ) */
            ctx.rule = rules.get(90);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("map", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Inputs", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_INPUT);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COLON);
            tree.add(next);
            subtree = parse__gen26(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "call_input",
            current,
            nonterminal_first.get(72),
            rules.get(90)
        ));
    }
    public ParseTree parse__gen26(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen26(ctx);
    }
    private static ParseTree parse__gen26(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[15][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(73, "_gen26"));
        ctx.nonterminal = "_gen26";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(73).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(73).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 86) {
            /* $_gen26 = $mapping $_gen27 */
            ctx.rule = rules.get(86);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_mapping(ctx);
            tree.add(subtree);
            subtree = parse__gen27(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_document(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_document(ctx);
    }
    private static ParseTree parse_document(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[16][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(74, "document"));
        ctx.nonterminal = "document";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(74).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(74).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 8) {
            /* $document = $_gen0 $_gen2 -> Document( imports=$0, definitions=$1 ) */
            ctx.rule = rules.get(8);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("imports", 0);
            parameters.put("definitions", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("Document", parameters));
            subtree = parse__gen0(ctx);
            tree.add(subtree);
            subtree = parse__gen2(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen34(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen34(ctx);
    }
    private static ParseTree parse__gen34(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[17][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(75, "_gen34"));
        ctx.nonterminal = "_gen34";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(75).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(75).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 131) {
            /* $_gen34 = :comma $object_kv $_gen34 */
            ctx.rule = rules.get(131);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COMMA);
            tree.add(next);
            tree.setListSeparator(next);
            subtree = parse_object_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen34(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_cmd_param_kv(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_cmd_param_kv(ctx);
    }
    private static ParseTree parse_cmd_param_kv(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[18][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(76, "cmd_param_kv"));
        ctx.nonterminal = "cmd_param_kv";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "cmd_param_kv",
                nonterminal_first.get(76),
                nonterminal_rules.get(76)
            ));
        }
        if (rule == 47) {
            /* $cmd_param_kv = :cmd_attr_hint :identifier :equal $e -> CommandParameterAttr( key=$1, value=$3 ) */
            ctx.rule = rules.get(47);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("key", 1);
            parameters.put("value", 3);
            tree.setAstTransformation(new AstTransformNodeCreator("CommandParameterAttr", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_EQUAL);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "cmd_param_kv",
            current,
            nonterminal_first.get(76),
            rules.get(47)
        ));
    }
    public ParseTree parse__gen23(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen23(ctx);
    }
    private static ParseTree parse__gen23(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[20][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(78, "_gen23"));
        ctx.nonterminal = "_gen23";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(78).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(78).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 78) {
            /* $_gen23 = $call_body */
            ctx.rule = rules.get(78);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_call_body(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_map(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_map(ctx);
    }
    private static ParseTree parse_map(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[21][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(79, "map"));
        ctx.nonterminal = "map";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "map",
                nonterminal_first.get(79),
                nonterminal_rules.get(79)
            ));
        }
        if (rule == 64) {
            /* $map = :lbrace $_gen18 :rbrace -> $1 */
            ctx.rule = rules.get(64);
            tree.setAstTransformation(new AstTransformSubstitution(1));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen18(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "map",
            current,
            nonterminal_first.get(79),
            rules.get(64)
        ));
    }
    public ParseTree parse_scatter(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_scatter(ctx);
    }
    private static ParseTree parse_scatter(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[22][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(80, "scatter"));
        ctx.nonterminal = "scatter";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "scatter",
                nonterminal_first.get(80),
                nonterminal_rules.get(80)
            ));
        }
        if (rule == 95) {
            /* $scatter = :scatter :lparen :identifier :in $e :rparen :lbrace $_gen20 :rbrace -> Scatter( item=$2, collection=$4, body=$7 ) */
            ctx.rule = rules.get(95);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("item", 2);
            parameters.put("collection", 4);
            parameters.put("body", 7);
            tree.setAstTransformation(new AstTransformNodeCreator("Scatter", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_SCATTER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LPAREN);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IN);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RPAREN);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen20(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "scatter",
            current,
            nonterminal_first.get(80),
            rules.get(95)
        ));
    }
    public ParseTree parse_sections(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_sections(ctx);
    }
    private static ParseTree parse_sections(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[23][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(81, "sections"));
        ctx.nonterminal = "sections";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "sections",
                nonterminal_first.get(81),
                nonterminal_rules.get(81)
            ));
        }
        if (rule == 24) {
            /* $sections = $command */
            ctx.rule = rules.get(24);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_command(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 25) {
            /* $sections = $outputs */
            ctx.rule = rules.get(25);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_outputs(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 26) {
            /* $sections = $runtime */
            ctx.rule = rules.get(26);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_runtime(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 27) {
            /* $sections = $parameter_meta */
            ctx.rule = rules.get(27);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_parameter_meta(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 28) {
            /* $sections = $meta */
            ctx.rule = rules.get(28);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_meta(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "sections",
            current,
            nonterminal_first.get(81),
            rules.get(28)
        ));
    }
    public ParseTree parse__gen0(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen0(ctx);
    }
    private static ParseTree parse__gen0(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[25][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(83, "_gen0"));
        ctx.nonterminal = "_gen0";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(83).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(83).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 0) {
            /* $_gen0 = $import $_gen1 */
            ctx.rule = rules.get(0);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_import(ctx);
            tree.add(subtree);
            subtree = parse__gen1(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen17(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen17(ctx);
    }
    private static ParseTree parse__gen17(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[26][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(84, "_gen17"));
        ctx.nonterminal = "_gen17";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(84).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(84).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 52) {
            /* $_gen17 = $output_kv $_gen17 */
            ctx.rule = rules.get(52);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_output_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen17(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_meta(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_meta(ctx);
    }
    private static ParseTree parse_meta(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[27][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(85, "meta"));
        ctx.nonterminal = "meta";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "meta",
                nonterminal_first.get(85),
                nonterminal_rules.get(85)
            ));
        }
        if (rule == 59) {
            /* $meta = :meta $map -> Meta( map=$1 ) */
            ctx.rule = rules.get(59);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("map", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("Meta", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_META);
            tree.add(next);
            subtree = parse_map(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "meta",
            current,
            nonterminal_first.get(85),
            rules.get(59)
        ));
    }
    public ParseTree parse__gen27(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen27(ctx);
    }
    private static ParseTree parse__gen27(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[28][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(86, "_gen27"));
        ctx.nonterminal = "_gen27";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(86).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(86).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 87) {
            /* $_gen27 = :comma $mapping $_gen27 */
            ctx.rule = rules.get(87);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COMMA);
            tree.add(next);
            tree.setListSeparator(next);
            subtree = parse_mapping(ctx);
            tree.add(subtree);
            subtree = parse__gen27(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen30(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen30(ctx);
    }
    private static ParseTree parse__gen30(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[29][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(87, "_gen30"));
        ctx.nonterminal = "_gen30";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(87).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(87).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 102) {
            /* $_gen30 = :comma $type_e $_gen30 */
            ctx.rule = rules.get(102);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COMMA);
            tree.add(next);
            tree.setListSeparator(next);
            subtree = parse_type_e(ctx);
            tree.add(subtree);
            subtree = parse__gen30(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen11(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen11(ctx);
    }
    private static ParseTree parse__gen11(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[30][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(88, "_gen11"));
        ctx.nonterminal = "_gen11";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(88).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(88).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 36) {
            /* $_gen11 = $cmd_param_kv $_gen12 */
            ctx.rule = rules.get(36);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_cmd_param_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen12(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_wf_body_element(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_wf_body_element(ctx);
    }
    private static ParseTree parse_wf_body_element(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[31][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(89, "wf_body_element"));
        ctx.nonterminal = "wf_body_element";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "wf_body_element",
                nonterminal_first.get(89),
                nonterminal_rules.get(89)
            ));
        }
        if (rule == 71) {
            /* $wf_body_element = $call */
            ctx.rule = rules.get(71);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_call(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 72) {
            /* $wf_body_element = $declaration */
            ctx.rule = rules.get(72);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_declaration(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 73) {
            /* $wf_body_element = $while_loop */
            ctx.rule = rules.get(73);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_while_loop(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 74) {
            /* $wf_body_element = $if_stmt */
            ctx.rule = rules.get(74);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_if_stmt(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 75) {
            /* $wf_body_element = $scatter */
            ctx.rule = rules.get(75);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_scatter(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "wf_body_element",
            current,
            nonterminal_first.get(89),
            rules.get(75)
        ));
    }
    public ParseTree parse__gen12(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen12(ctx);
    }
    private static ParseTree parse__gen12(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[32][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(90, "_gen12"));
        ctx.nonterminal = "_gen12";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(90).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(90).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 37) {
            /* $_gen12 = $cmd_param_kv $_gen12 */
            ctx.rule = rules.get(37);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_cmd_param_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen12(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen21(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen21(ctx);
    }
    private static ParseTree parse__gen21(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[33][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(91, "_gen21"));
        ctx.nonterminal = "_gen21";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(91).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(91).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 67) {
            /* $_gen21 = $wf_body_element $_gen21 */
            ctx.rule = rules.get(67);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_wf_body_element(ctx);
            tree.add(subtree);
            subtree = parse__gen21(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_declaration(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_declaration(ctx);
    }
    private static ParseTree parse_declaration(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[34][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(92, "declaration"));
        ctx.nonterminal = "declaration";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "declaration",
                nonterminal_first.get(92),
                nonterminal_rules.get(92)
            ));
        }
        if (rule == 98) {
            /* $declaration = $type_e :identifier $_gen28 -> Declaration( type=$0, name=$1, expression=$2 ) */
            ctx.rule = rules.get(98);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("type", 0);
            parameters.put("name", 1);
            parameters.put("expression", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Declaration", parameters));
            subtree = parse_type_e(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            subtree = parse__gen28(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "declaration",
            current,
            nonterminal_first.get(92),
            rules.get(98)
        ));
    }
    public ParseTree parse__gen8(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen8(ctx);
    }
    private static ParseTree parse__gen8(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[35][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(93, "_gen8"));
        ctx.nonterminal = "_gen8";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(93).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(93).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 20) {
            /* $_gen8 = $sections $_gen8 */
            ctx.rule = rules.get(20);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_sections(ctx);
            tree.add(subtree);
            subtree = parse__gen8(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen4(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen4(ctx);
    }
    private static ParseTree parse__gen4(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[36][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(94, "_gen4"));
        ctx.nonterminal = "_gen4";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(94).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(94).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 11) {
            /* $_gen4 = $import_namespace */
            ctx.rule = rules.get(11);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_import_namespace(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen2(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen2(ctx);
    }
    private static ParseTree parse__gen2(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[37][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(95, "_gen2"));
        ctx.nonterminal = "_gen2";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(95).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(95).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 4) {
            /* $_gen2 = $workflow_or_task $_gen3 */
            ctx.rule = rules.get(4);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_workflow_or_task(ctx);
            tree.add(subtree);
            subtree = parse__gen3(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen22(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen22(ctx);
    }
    private static ParseTree parse__gen22(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[38][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(96, "_gen22"));
        ctx.nonterminal = "_gen22";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(96).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(96).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 76) {
            /* $_gen22 = $alias */
            ctx.rule = rules.get(76);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_alias(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_command_part(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_command_part(ctx);
    }
    private static ParseTree parse_command_part(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[39][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(97, "command_part"));
        ctx.nonterminal = "command_part";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "command_part",
                nonterminal_first.get(97),
                nonterminal_rules.get(97)
            ));
        }
        if (rule == 34) {
            /* $command_part = :cmd_part */
            ctx.rule = rules.get(34);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_CMD_PART);
            tree.add(next);
            return tree;
        }
        else if (rule == 35) {
            /* $command_part = $cmd_param */
            ctx.rule = rules.get(35);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_cmd_param(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "command_part",
            current,
            nonterminal_first.get(97),
            rules.get(35)
        ));
    }
    public ParseTree parse__gen6(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen6(ctx);
    }
    private static ParseTree parse__gen6(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[40][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(98, "_gen6"));
        ctx.nonterminal = "_gen6";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(98).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(98).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 16) {
            /* $_gen6 = $declaration $_gen6 */
            ctx.rule = rules.get(16);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_declaration(ctx);
            tree.add(subtree);
            subtree = parse__gen6(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen28(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen28(ctx);
    }
    private static ParseTree parse__gen28(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[41][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(99, "_gen28"));
        ctx.nonterminal = "_gen28";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(99).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(99).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 96) {
            /* $_gen28 = $setter */
            ctx.rule = rules.get(96);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_setter(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen3(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen3(ctx);
    }
    private static ParseTree parse__gen3(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[42][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(100, "_gen3"));
        ctx.nonterminal = "_gen3";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(100).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(100).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 5) {
            /* $_gen3 = $workflow_or_task $_gen3 */
            ctx.rule = rules.get(5);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_workflow_or_task(ctx);
            tree.add(subtree);
            subtree = parse__gen3(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_runtime(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_runtime(ctx);
    }
    private static ParseTree parse_runtime(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[43][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(101, "runtime"));
        ctx.nonterminal = "runtime";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "runtime",
                nonterminal_first.get(101),
                nonterminal_rules.get(101)
            ));
        }
        if (rule == 57) {
            /* $runtime = :runtime $map -> Runtime( map=$1 ) */
            ctx.rule = rules.get(57);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("map", 1);
            tree.setAstTransformation(new AstTransformNodeCreator("Runtime", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RUNTIME);
            tree.add(next);
            subtree = parse_map(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "runtime",
            current,
            nonterminal_first.get(101),
            rules.get(57)
        ));
    }
    public ParseTree parse__gen29(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen29(ctx);
    }
    private static ParseTree parse__gen29(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[44][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(102, "_gen29"));
        ctx.nonterminal = "_gen29";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(102).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(102).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 101) {
            /* $_gen29 = $type_e $_gen30 */
            ctx.rule = rules.get(101);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_type_e(ctx);
            tree.add(subtree);
            subtree = parse__gen30(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen33(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen33(ctx);
    }
    private static ParseTree parse__gen33(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[45][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(103, "_gen33"));
        ctx.nonterminal = "_gen33";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(103).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(103).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 130) {
            /* $_gen33 = $object_kv $_gen34 */
            ctx.rule = rules.get(130);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_object_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen34(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen9(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen9(ctx);
    }
    private static ParseTree parse__gen9(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[46][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(104, "_gen9"));
        ctx.nonterminal = "_gen9";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(104).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(104).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 29) {
            /* $_gen9 = $command_part $_gen10 */
            ctx.rule = rules.get(29);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_command_part(ctx);
            tree.add(subtree);
            subtree = parse__gen10(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_postfix_quantifier(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_postfix_quantifier(ctx);
    }
    private static ParseTree parse_postfix_quantifier(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[47][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(105, "postfix_quantifier"));
        ctx.nonterminal = "postfix_quantifier";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "postfix_quantifier",
                nonterminal_first.get(105),
                nonterminal_rules.get(105)
            ));
        }
        if (rule == 48) {
            /* $postfix_quantifier = :qmark */
            ctx.rule = rules.get(48);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_QMARK);
            tree.add(next);
            return tree;
        }
        else if (rule == 49) {
            /* $postfix_quantifier = :plus */
            ctx.rule = rules.get(49);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_PLUS);
            tree.add(next);
            return tree;
        }
        else if (rule == 50) {
            /* $postfix_quantifier = :asterisk */
            ctx.rule = rules.get(50);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_ASTERISK);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "postfix_quantifier",
            current,
            nonterminal_first.get(105),
            rules.get(50)
        ));
    }
    public ParseTree parse__gen24(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen24(ctx);
    }
    private static ParseTree parse__gen24(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[48][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(106, "_gen24"));
        ctx.nonterminal = "_gen24";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(106).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(106).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 81) {
            /* $_gen24 = $call_input $_gen25 */
            ctx.rule = rules.get(81);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_call_input(ctx);
            tree.add(subtree);
            subtree = parse__gen25(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen32(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen32(ctx);
    }
    private static ParseTree parse__gen32(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[49][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(107, "_gen32"));
        ctx.nonterminal = "_gen32";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(107).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(107).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 124) {
            /* $_gen32 = :comma $e $_gen32 */
            ctx.rule = rules.get(124);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COMMA);
            tree.add(next);
            tree.setListSeparator(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            subtree = parse__gen32(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_task(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_task(ctx);
    }
    private static ParseTree parse_task(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[50][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(108, "task"));
        ctx.nonterminal = "task";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "task",
                nonterminal_first.get(108),
                nonterminal_rules.get(108)
            ));
        }
        if (rule == 23) {
            /* $task = :task :identifier :lbrace $_gen5 $_gen7 :rbrace -> Task( name=$1, declarations=$3, sections=$4 ) */
            ctx.rule = rules.get(23);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("name", 1);
            parameters.put("declarations", 3);
            parameters.put("sections", 4);
            tree.setAstTransformation(new AstTransformNodeCreator("Task", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_TASK);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen5(ctx);
            tree.add(subtree);
            subtree = parse__gen7(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "task",
            current,
            nonterminal_first.get(108),
            rules.get(23)
        ));
    }
    public ParseTree parse_outputs(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_outputs(ctx);
    }
    private static ParseTree parse_outputs(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[51][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(109, "outputs"));
        ctx.nonterminal = "outputs";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "outputs",
                nonterminal_first.get(109),
                nonterminal_rules.get(109)
            ));
        }
        if (rule == 55) {
            /* $outputs = :output :lbrace $_gen16 :rbrace -> Outputs( attributes=$2 ) */
            ctx.rule = rules.get(55);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("attributes", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("Outputs", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_OUTPUT);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen16(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "outputs",
            current,
            nonterminal_first.get(109),
            rules.get(55)
        ));
    }
    public ParseTree parse_workflow(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_workflow(ctx);
    }
    private static ParseTree parse_workflow(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[52][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(110, "workflow"));
        ctx.nonterminal = "workflow";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "workflow",
                nonterminal_first.get(110),
                nonterminal_rules.get(110)
            ));
        }
        if (rule == 70) {
            /* $workflow = :workflow :identifier :lbrace $_gen20 :rbrace -> Workflow( name=$1, body=$3 ) */
            ctx.rule = rules.get(70);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("name", 1);
            parameters.put("body", 3);
            tree.setAstTransformation(new AstTransformNodeCreator("Workflow", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_WORKFLOW);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen20(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "workflow",
            current,
            nonterminal_first.get(110),
            rules.get(70)
        ));
    }
    public ParseTree parse_import_namespace(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_import_namespace(ctx);
    }
    private static ParseTree parse_import_namespace(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[53][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(111, "import_namespace"));
        ctx.nonterminal = "import_namespace";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "import_namespace",
                nonterminal_first.get(111),
                nonterminal_rules.get(111)
            ));
        }
        if (rule == 14) {
            /* $import_namespace = :as :identifier -> $1 */
            ctx.rule = rules.get(14);
            tree.setAstTransformation(new AstTransformSubstitution(1));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_AS);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "import_namespace",
            current,
            nonterminal_first.get(111),
            rules.get(14)
        ));
    }
    public ParseTree parse_output_kv(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_output_kv(ctx);
    }
    private static ParseTree parse_output_kv(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[54][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(112, "output_kv"));
        ctx.nonterminal = "output_kv";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "output_kv",
                nonterminal_first.get(112),
                nonterminal_rules.get(112)
            ));
        }
        if (rule == 56) {
            /* $output_kv = $type_e :identifier :equal $e -> Output( type=$0, var=$1, expression=$3 ) */
            ctx.rule = rules.get(56);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("type", 0);
            parameters.put("var", 1);
            parameters.put("expression", 3);
            tree.setAstTransformation(new AstTransformNodeCreator("Output", parameters));
            subtree = parse_type_e(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_EQUAL);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "output_kv",
            current,
            nonterminal_first.get(112),
            rules.get(56)
        ));
    }
    public ParseTree parse_call_body(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_call_body(ctx);
    }
    private static ParseTree parse_call_body(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[55][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(113, "call_body"));
        ctx.nonterminal = "call_body";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "call_body",
                nonterminal_first.get(113),
                nonterminal_rules.get(113)
            ));
        }
        if (rule == 85) {
            /* $call_body = :lbrace $_gen5 $_gen24 :rbrace -> CallBody( declarations=$1, io=$2 ) */
            ctx.rule = rules.get(85);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("declarations", 1);
            parameters.put("io", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("CallBody", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen5(ctx);
            tree.add(subtree);
            subtree = parse__gen24(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "call_body",
            current,
            nonterminal_first.get(113),
            rules.get(85)
        ));
    }
    public ParseTree parse_mapping(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_mapping(ctx);
    }
    private static ParseTree parse_mapping(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[56][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(114, "mapping"));
        ctx.nonterminal = "mapping";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "mapping",
                nonterminal_first.get(114),
                nonterminal_rules.get(114)
            ));
        }
        if (rule == 91) {
            /* $mapping = :identifier :equal $e -> IOMapping( key=$0, value=$2 ) */
            ctx.rule = rules.get(91);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("key", 0);
            parameters.put("value", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("IOMapping", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_EQUAL);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "mapping",
            current,
            nonterminal_first.get(114),
            rules.get(91)
        ));
    }
    public ParseTree parse_alias(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_alias(ctx);
    }
    private static ParseTree parse_alias(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[57][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(115, "alias"));
        ctx.nonterminal = "alias";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "alias",
                nonterminal_first.get(115),
                nonterminal_rules.get(115)
            ));
        }
        if (rule == 92) {
            /* $alias = :as :identifier -> $1 */
            ctx.rule = rules.get(92);
            tree.setAstTransformation(new AstTransformSubstitution(1));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_AS);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "alias",
            current,
            nonterminal_first.get(115),
            rules.get(92)
        ));
    }
    public ParseTree parse_kv(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_kv(ctx);
    }
    private static ParseTree parse_kv(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[58][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(116, "kv"));
        ctx.nonterminal = "kv";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "kv",
                nonterminal_first.get(116),
                nonterminal_rules.get(116)
            ));
        }
        if (rule == 65) {
            /* $kv = :identifier :colon $e -> RuntimeAttribute( key=$0, value=$2 ) */
            ctx.rule = rules.get(65);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("key", 0);
            parameters.put("value", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("RuntimeAttribute", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COLON);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "kv",
            current,
            nonterminal_first.get(116),
            rules.get(65)
        ));
    }
    public ParseTree parse_command(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_command(ctx);
    }
    private static ParseTree parse_command(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[59][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(117, "command"));
        ctx.nonterminal = "command";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "command",
                nonterminal_first.get(117),
                nonterminal_rules.get(117)
            ));
        }
        if (rule == 33) {
            /* $command = :raw_command :raw_cmd_start $_gen9 :raw_cmd_end -> RawCommand( parts=$2 ) */
            ctx.rule = rules.get(33);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("parts", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("RawCommand", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RAW_COMMAND);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RAW_CMD_START);
            tree.add(next);
            subtree = parse__gen9(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RAW_CMD_END);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "command",
            current,
            nonterminal_first.get(117),
            rules.get(33)
        ));
    }
    public ParseTree parse_while_loop(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_while_loop(ctx);
    }
    private static ParseTree parse_while_loop(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[60][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(118, "while_loop"));
        ctx.nonterminal = "while_loop";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "while_loop",
                nonterminal_first.get(118),
                nonterminal_rules.get(118)
            ));
        }
        if (rule == 93) {
            /* $while_loop = :while :lparen $e :rparen :lbrace $_gen20 :rbrace -> WhileLoop( expression=$2, body=$5 ) */
            ctx.rule = rules.get(93);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("expression", 2);
            parameters.put("body", 5);
            tree.setAstTransformation(new AstTransformNodeCreator("WhileLoop", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_WHILE);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LPAREN);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RPAREN);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen20(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "while_loop",
            current,
            nonterminal_first.get(118),
            rules.get(93)
        ));
    }
    public ParseTree parse_workflow_or_task(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_workflow_or_task(ctx);
    }
    private static ParseTree parse_workflow_or_task(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[61][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(119, "workflow_or_task"));
        ctx.nonterminal = "workflow_or_task";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "workflow_or_task",
                nonterminal_first.get(119),
                nonterminal_rules.get(119)
            ));
        }
        if (rule == 9) {
            /* $workflow_or_task = $workflow */
            ctx.rule = rules.get(9);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_workflow(ctx);
            tree.add(subtree);
            return tree;
        }
        else if (rule == 10) {
            /* $workflow_or_task = $task */
            ctx.rule = rules.get(10);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_task(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "workflow_or_task",
            current,
            nonterminal_first.get(119),
            rules.get(10)
        ));
    }
    public ParseTree parse__gen16(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen16(ctx);
    }
    private static ParseTree parse__gen16(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[62][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(120, "_gen16"));
        ctx.nonterminal = "_gen16";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(120).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(120).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 51) {
            /* $_gen16 = $output_kv $_gen17 */
            ctx.rule = rules.get(51);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_output_kv(ctx);
            tree.add(subtree);
            subtree = parse__gen17(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse_if_stmt(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_if_stmt(ctx);
    }
    private static ParseTree parse_if_stmt(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[63][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(121, "if_stmt"));
        ctx.nonterminal = "if_stmt";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "if_stmt",
                nonterminal_first.get(121),
                nonterminal_rules.get(121)
            ));
        }
        if (rule == 94) {
            /* $if_stmt = :if :lparen $e :rparen :lbrace $_gen20 :rbrace -> If( expression=$2, body=$5 ) */
            ctx.rule = rules.get(94);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("expression", 2);
            parameters.put("body", 5);
            tree.setAstTransformation(new AstTransformNodeCreator("If", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IF);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LPAREN);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RPAREN);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_LBRACE);
            tree.add(next);
            subtree = parse__gen20(ctx);
            tree.add(subtree);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_RBRACE);
            tree.add(next);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "if_stmt",
            current,
            nonterminal_first.get(121),
            rules.get(94)
        ));
    }
    public ParseTree parse_object_kv(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse_object_kv(ctx);
    }
    private static ParseTree parse_object_kv(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[64][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(122, "object_kv"));
        ctx.nonterminal = "object_kv";
        tree.setList(false);
        if (current == null) {
            throw new SyntaxError(ctx.error_formatter.unexpectedEof(
                "object_kv",
                nonterminal_first.get(122),
                nonterminal_rules.get(122)
            ));
        }
        if (rule == 100) {
            /* $object_kv = :identifier :colon $e -> ObjectKV( key=$0, value=$2 ) */
            ctx.rule = rules.get(100);
            LinkedHashMap<String, Integer> parameters = new LinkedHashMap<String, Integer>();
            parameters.put("key", 0);
            parameters.put("value", 2);
            tree.setAstTransformation(new AstTransformNodeCreator("ObjectKV", parameters));
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_IDENTIFIER);
            tree.add(next);
            next = expect(ctx, WdlTerminalIdentifier.TERMINAL_COLON);
            tree.add(next);
            subtree = parse_e(ctx);
            tree.add(subtree);
            return tree;
        }
        throw new SyntaxError(ctx.error_formatter.unexpectedSymbol(
            "object_kv",
            current,
            nonterminal_first.get(122),
            rules.get(100)
        ));
    }
    public ParseTree parse__gen31(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen31(ctx);
    }
    private static ParseTree parse__gen31(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[65][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(123, "_gen31"));
        ctx.nonterminal = "_gen31";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(123).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(123).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 123) {
            /* $_gen31 = $e $_gen32 */
            ctx.rule = rules.get(123);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_e(ctx);
            tree.add(subtree);
            subtree = parse__gen32(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen10(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen10(ctx);
    }
    private static ParseTree parse__gen10(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[66][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(124, "_gen10"));
        ctx.nonterminal = "_gen10";
        tree.setList(true);
        if ( current != null &&
             !nonterminal_first.get(124).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(124).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 30) {
            /* $_gen10 = $command_part $_gen10 */
            ctx.rule = rules.get(30);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_command_part(ctx);
            tree.add(subtree);
            subtree = parse__gen10(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    public ParseTree parse__gen15(List<Terminal> tokens, SyntaxErrorFormatter error_formatter) throws SyntaxError {
        ParserContext ctx = new ParserContext(new TokenStream(tokens), error_formatter);
        return parse__gen15(ctx);
    }
    private static ParseTree parse__gen15(ParserContext ctx) throws SyntaxError {
        Terminal current = ctx.tokens.current();
        Terminal next;
        ParseTree subtree;
        int rule = (current != null) ? table[67][current.getId()] : -1;
        ParseTree tree = new ParseTree( new NonTerminal(125, "_gen15"));
        ctx.nonterminal = "_gen15";
        tree.setList(false);
        if ( current != null &&
             !nonterminal_first.get(125).contains(terminal_map.get(current.getId())) &&
              nonterminal_follow.get(125).contains(terminal_map.get(current.getId())) ) {
            return tree;
        }
        if (current == null) {
            return tree;
        }
        if (rule == 44) {
            /* $_gen15 = $postfix_quantifier */
            ctx.rule = rules.get(44);
            tree.setAstTransformation(new AstTransformSubstitution(0));
            subtree = parse_postfix_quantifier(ctx);
            tree.add(subtree);
            return tree;
        }
        return tree;
    }
    /* Section: Lexer */
    private Map<String, List<HermesRegex>> regex = null;
    private interface LexerOutput {}
    private class LexerRegexOutput implements LexerOutput {
        public WdlTerminalIdentifier terminal;
        public int group;
        public Method function;
        LexerRegexOutput(WdlTerminalIdentifier terminal, int group, Method function) {
            this.terminal = terminal;
            this.group = group;
            this.function = function;
        }
        public String toString() {
            return String.format("<LexerRegexOutput terminal=%s, group=%d, func=%s>", this.terminal, this.group, this.function);
        }
    }
    private class LexerStackPush implements LexerOutput {
        public String mode;
        LexerStackPush(String mode) {
            this.mode = mode;
        }
    }
    private class LexerAction implements LexerOutput {
        public String action;
        LexerAction(String action) {
            this.action = action;
        }
    }
    private class HermesRegex {
        public Pattern pattern;
        public List<LexerOutput> outputs;
        HermesRegex(Pattern pattern, List<LexerOutput> outputs) {
            this.pattern = pattern;
            this.outputs = outputs;
        }
        public String toString() {
            return String.format("<HermesRegex pattern=%s, outputs=%s>", this.pattern, this.outputs);
        }
    }
    private class LineColumn {
        public int line, col;
        public LineColumn(int line, int col) {
            this.line = line;
            this.col = col;
        }
        public String toString() {
            return String.format("<LineColumn: line=%d column=%d>", this.line, this.col);
        }
    }
    private class LexerContext {
        public String string;
        public String resource;
        public int line;
        public int col;
        public Stack<String> stack;
        public Object context;
        public List<Terminal> terminals;
        LexerContext(String string, String resource) {
            this.string = string;
            this.resource = resource;
            this.line = 1;
            this.col = 1;
            this.stack = new Stack<String>();
            this.stack.push("default");
            this.terminals = new ArrayList<Terminal>();
        }
        public void advance(String match) {
            LineColumn lc = advance_line_col(match, match.length());
            this.line = lc.line;
            this.col = lc.col;
            this.string = this.string.substring(match.length());
        }
        public LineColumn advance_line_col(String match, int length) {
            LineColumn lc = new LineColumn(this.line, this.col);
            for (int i = 0; i < length && i < match.length(); i++) {
                if (match.charAt(i) == '\n') {
                    lc.line += 1;
                    lc.col = 1;
                } else {
                    lc.col += 1;
                }
            }
            return lc;
        }
    }
    private void emit(LexerContext lctx, TerminalIdentifier terminal, String source_string, int line, int col) {
        lctx.terminals.add(new Terminal(terminal.id(), terminal.string(), source_string, lctx.resource, line, col));
    }
    /**
     * The default function that is called on every regex match during lexical analysis.
     * By default, this simply calls the emit() function with all of the same parameters.
     * This can be overridden in the grammar file to provide a different default action.
     *
     * @param lctx The current state of the lexical analyzer
     * @param terminal The current terminal that was matched
     * @param source_string The source code that was matched
     * @param line The line where the match happened
     * @param col The column where the match happened
     * @return void
     */
    public void default_action(LexerContext lctx, TerminalIdentifier terminal, String source_string, int line, int col) {
        emit(lctx, terminal, source_string, line, col);
    }
    /* START USER CODE */
    /* END USER CODE */
    public Object init() {
        return null;
    }
    public void destroy(Object context) {
        return;
    }
    private Method getFunction(String name) throws SyntaxError {
        try {
            return getClass().getMethod(
                name,
                LexerContext.class,
                TerminalIdentifier.class,
                String.class,
                int.class,
                int.class
            );
        } catch (NoSuchMethodException e) {
            throw new SyntaxError("No such method: " + name);
        }
    }
    private void lexer_init() throws SyntaxError {
        this.regex = new HashMap<String, List<HermesRegex>>();
        this.regex.put("default", Arrays.asList(new HermesRegex[] {
            new HermesRegex(
                Pattern.compile("\\s+"),
                Arrays.asList(new LexerOutput[] {
                })
            ),
            new HermesRegex(
                Pattern.compile("/\\*(.*?)\\*/", Pattern.DOTALL),
                Arrays.asList(new LexerOutput[] {
                })
            ),
            new HermesRegex(
                Pattern.compile("#.*"),
                Arrays.asList(new LexerOutput[] {
                })
            ),
            new HermesRegex(
                Pattern.compile("task(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_TASK,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("(call)\\s+"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CALL,
                        1,
                        getFunction("default_action")
                    ),
                    new LexerStackPush("task_fqn"),
                })
            ),
            new HermesRegex(
                Pattern.compile("workflow(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_WORKFLOW,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("import(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IMPORT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("input(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_INPUT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("output(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_OUTPUT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("as(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_AS,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("if(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IF,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("while(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_WHILE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("runtime(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RUNTIME,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("scatter(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_SCATTER,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerStackPush("scatter"),
                })
            ),
            new HermesRegex(
                Pattern.compile("command\\s*(?=<<<)"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerStackPush("raw_command2"),
                })
            ),
            new HermesRegex(
                Pattern.compile("command\\s*(?=\\{)"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RAW_COMMAND,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerStackPush("raw_command"),
                })
            ),
            new HermesRegex(
                Pattern.compile("parameter_meta(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_PARAMETER_META,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("meta(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_META,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("(true|false)(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_BOOLEAN,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("(object)\\s*(\\{)"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_OBJECT,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LBRACE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("(Array|Map|Object|Boolean|Int|Float|Uri|File|String)(?![a-zA-Z0-9_])(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_TYPE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("[a-zA-Z]([a-zA-Z0-9_])*"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\"([^\"]+)\""),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_STRING,
                        1,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("'([^']+)'"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_STRING,
                        1,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile(":"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_COLON,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile(","),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_COMMA,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("=="),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_DOUBLE_EQUAL,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\|\\|"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_DOUBLE_PIPE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\&\\&"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_DOUBLE_AMPERSAND,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("!="),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_NOT_EQUAL,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("="),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_EQUAL,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\."),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_DOT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\{"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LBRACE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\}"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RBRACE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\("),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LPAREN,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\)"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RPAREN,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\["),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LSQUARE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\]"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RSQUARE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\+"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_PLUS,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\*"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_ASTERISK,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("-"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_DASH,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("/"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_SLASH,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("%"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_PERCENT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("<="),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LTEQ,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("<"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile(">="),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_GTEQ,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile(">"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_GT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("!"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_NOT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("-?[0-9]+\\.[0-9]+"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_FLOAT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("[0-9]+"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_INTEGER,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
        }));
        this.regex.put("task_fqn", Arrays.asList(new HermesRegex[] {
            new HermesRegex(
                Pattern.compile("\\s+"),
                Arrays.asList(new LexerOutput[] {
                })
            ),
            new HermesRegex(
                Pattern.compile("[a-zA-Z]([a-zA-Z0-9_])*(\\.[a-zA-Z]([a-zA-Z0-9_])*)*"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_NS_IDENTIFIER,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerAction("pop"),
                })
            ),
        }));
        this.regex.put("scatter", Arrays.asList(new HermesRegex[] {
            new HermesRegex(
                Pattern.compile("\\s+"),
                Arrays.asList(new LexerOutput[] {
                })
            ),
            new HermesRegex(
                Pattern.compile("\\)"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RPAREN,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerAction("pop"),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\("),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LPAREN,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\."),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_DOT,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\["),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LSQUARE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\]"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RSQUARE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("in(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IN,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("[a-zA-Z]([a-zA-Z0-9_])*"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
        }));
        this.regex.put("raw_command", Arrays.asList(new HermesRegex[] {
            new HermesRegex(
                Pattern.compile("\\{"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RAW_CMD_START,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\}"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RAW_CMD_END,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerAction("pop"),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\$\\{"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerStackPush("cmd_param"),
                })
            ),
            new HermesRegex(
                Pattern.compile("(.*?)(?=\\$\\{|\\})", Pattern.DOTALL),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CMD_PART,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
        }));
        this.regex.put("raw_command2", Arrays.asList(new HermesRegex[] {
            new HermesRegex(
                Pattern.compile("<<<"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RAW_CMD_START,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile(">>>"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RAW_CMD_END,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerAction("pop"),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\$\\{"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CMD_PARAM_START,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerStackPush("cmd_param"),
                })
            ),
            new HermesRegex(
                Pattern.compile("(.*?)(?=\\$\\{|>>>)", Pattern.DOTALL),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CMD_PART,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
        }));
        this.regex.put("cmd_param", Arrays.asList(new HermesRegex[] {
            new HermesRegex(
                Pattern.compile("\\s+"),
                Arrays.asList(new LexerOutput[] {
                })
            ),
            new HermesRegex(
                Pattern.compile("\\}"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CMD_PARAM_END,
                        0,
                        getFunction("default_action")
                    ),
                    new LexerAction("pop"),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\["),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_LSQUARE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\]"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_RSQUARE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("="),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_EQUAL,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\?"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_QMARK,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\+"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_PLUS,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\\*"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_ASTERISK,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("[0-9]+"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_INTEGER,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("(true|false)(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_BOOLEAN,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("(Array|Map|Object|Boolean|Int|Float|Uri|File|String)(?![a-zA-Z0-9_])(?![a-zA-Z0-9_])"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_TYPE,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("[a-zA-Z]([a-zA-Z0-9_])*(?=\\s*=)"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_CMD_ATTR_HINT,
                        -1,
                        getFunction("default_action")
                    ),
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("[a-zA-Z]([a-zA-Z0-9_])*"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_IDENTIFIER,
                        0,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("\"([^\"]+)\""),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_STRING,
                        1,
                        getFunction("default_action")
                    ),
                })
            ),
            new HermesRegex(
                Pattern.compile("'([^']+)'"),
                Arrays.asList(new LexerOutput[] {
                    new LexerRegexOutput(
                        WdlTerminalIdentifier.TERMINAL_STRING,
                        1,
                        getFunction("default_action")
                    ),
                })
            ),
        }));
    }
    private void unrecognized_token(String string, int line, int col) throws SyntaxError {
        String[] a = string.split("\n");
        String bad_line = string.split("\n")[line-1];
        StringBuffer spaces = new StringBuffer();
        for (int i = 0; i < col-1; i++) {
          spaces.append(' ');
        }
        String message = String.format(
            "Unrecognized token on line %d, column %d:\n\n%s\n%s^",
            line, col, bad_line, spaces
        );
        throw new SyntaxError(message);
    }
    private int next(LexerContext lctx) throws SyntaxError {
        String mode = lctx.stack.peek();
        for (int i = 0; i < this.regex.get(mode).size(); i++) {
            HermesRegex regex = this.regex.get(mode).get(i);
            Matcher matcher = regex.pattern.matcher(lctx.string);
            if (matcher.lookingAt()) {
                for (LexerOutput output : regex.outputs) {
                    if (output instanceof LexerStackPush) {
                        lctx.stack.push(((LexerStackPush) output).mode);
                    } else if (output instanceof LexerAction) {
                        LexerAction action = (LexerAction) output;
                        if (!action.action.equals("pop")) {
                            throw new SyntaxError("Invalid action");
                        }
                        if (action.action.equals("pop")) {
                            if (lctx.stack.empty()) {
                                throw new SyntaxError("Stack empty, cannot pop");
                            }
                            lctx.stack.pop();
                        }
                    } else if (output instanceof LexerRegexOutput) {
                        LexerRegexOutput regex_output = (LexerRegexOutput) output;
                        int group_line = lctx.line;
                        int group_col = lctx.col;
                        if (regex_output.group > 0) {
                            LineColumn lc = lctx.advance_line_col(matcher.group(0), matcher.start(regex_output.group));
                            group_line = lc.line;
                            group_col = lc.col;
                        }
                        try {
                            String source_string = (regex_output.group >= 0) ? matcher.group(regex_output.group) : "";
                            regex_output.function.invoke(
                                this,
                                lctx,
                                regex_output.terminal,
                                source_string,
                                group_line,
                                group_col
                            );
                        } catch (Exception e) {
                            e.printStackTrace();
                            throw new SyntaxError("Invalid method: " + regex_output.function);
                        }
                    }
                }
                lctx.advance(matcher.group(0));
                return matcher.group(0).length();
            }
        }
        return 0;
    }
    /**
     * Lexically analyze WDL source code, return a sequence of tokens.  Output of this
     * method should be used to construct a TerminalStream and then pass that to parse()
     *
     * @param string The WDL source code to analyze
     * @param resource A descriptor of where this code came from (usually a file path)
     * @return List of Terminal objects.
     * @throws SyntaxError If part of the source code could not lexically analyzed
     */
    public List<Terminal> lex(String string, String resource) throws SyntaxError {
        LexerContext lctx = new LexerContext(string, resource);
        Object context = this.init();
        String string_copy = new String(string);
        if (this.regex == null) {
            lexer_init();
        }
        while (lctx.string.length() > 0) {
            int match_length = this.next(lctx);
            if (match_length == 0) {
                this.unrecognized_token(string_copy, lctx.line, lctx.col);
            }
        }
        this.destroy(context);
        return lctx.terminals;
    }
    /* Section: Main */
}
