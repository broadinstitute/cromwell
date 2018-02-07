package wdl4s.parser;

public interface SharedWdlParserPieces {
    public static class SyntaxError extends Exception {
        public SyntaxError(String message) {
            super(message);
        }
    }
}
