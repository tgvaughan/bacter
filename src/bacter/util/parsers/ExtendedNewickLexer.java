// Generated from /home/tvaughan/code/beast_and_friends/bacter/src/bacter/util/parsers/ExtendedNewick.g4 by ANTLR 4.5.3
package bacter.util.parsers;
import org.antlr.v4.runtime.Lexer;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.Token;
import org.antlr.v4.runtime.TokenStream;
import org.antlr.v4.runtime.*;
import org.antlr.v4.runtime.atn.*;
import org.antlr.v4.runtime.dfa.DFA;
import org.antlr.v4.runtime.misc.*;

@SuppressWarnings({"all", "warnings", "unchecked", "unused", "cast"})
public class ExtendedNewickLexer extends Lexer {
	static { RuntimeMetaData.checkVersion("4.5.3", RuntimeMetaData.VERSION); }

	protected static final DFA[] _decisionToDFA;
	protected static final PredictionContextCache _sharedContextCache =
		new PredictionContextCache();
	public static final int
		T__0=1, T__1=2, T__2=3, T__3=4, T__4=5, T__5=6, T__6=7, T__7=8, T__8=9, 
		T__9=10, T__10=11, FLOAT=12, INT=13, RECTYPE=14, STRINGPRIM=15, WHITESPACE=16;
	public static String[] modeNames = {
		"DEFAULT_MODE"
	};

	public static final String[] ruleNames = {
		"T__0", "T__1", "T__2", "T__3", "T__4", "T__5", "T__6", "T__7", "T__8", 
		"T__9", "T__10", "FLOAT", "INT", "NNINT", "NZD", "D", "RECTYPE", "STRINGPRIM", 
		"WHITESPACE"
	};

	private static final String[] _LITERAL_NAMES = {
		null, "';'", "'('", "','", "')'", "':'", "'#'", "'[&'", "']'", "'='", 
		"'{'", "'}'"
	};
	private static final String[] _SYMBOLIC_NAMES = {
		null, null, null, null, null, null, null, null, null, null, null, null, 
		"FLOAT", "INT", "RECTYPE", "STRINGPRIM", "WHITESPACE"
	};
	public static final Vocabulary VOCABULARY = new VocabularyImpl(_LITERAL_NAMES, _SYMBOLIC_NAMES);

	/**
	 * @deprecated Use {@link #VOCABULARY} instead.
	 */
	@Deprecated
	public static final String[] tokenNames;
	static {
		tokenNames = new String[_SYMBOLIC_NAMES.length];
		for (int i = 0; i < tokenNames.length; i++) {
			tokenNames[i] = VOCABULARY.getLiteralName(i);
			if (tokenNames[i] == null) {
				tokenNames[i] = VOCABULARY.getSymbolicName(i);
			}

			if (tokenNames[i] == null) {
				tokenNames[i] = "<INVALID>";
			}
		}
	}

	@Override
	@Deprecated
	public String[] getTokenNames() {
		return tokenNames;
	}

	@Override

	public Vocabulary getVocabulary() {
		return VOCABULARY;
	}


	public ExtendedNewickLexer(CharStream input) {
		super(input);
		_interp = new LexerATNSimulator(this,_ATN,_decisionToDFA,_sharedContextCache);
	}

	@Override
	public String getGrammarFileName() { return "ExtendedNewick.g4"; }

	@Override
	public String[] getRuleNames() { return ruleNames; }

	@Override
	public String getSerializedATN() { return _serializedATN; }

	@Override
	public String[] getModeNames() { return modeNames; }

	@Override
	public ATN getATN() { return _ATN; }

	public static final String _serializedATN =
		"\3\u0430\ud6d1\u8206\uad2d\u4417\uaef1\u8d80\uaadd\2\22\u008d\b\1\4\2"+
		"\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7\4\b\t\b\4\t\t\t\4\n\t\n\4"+
		"\13\t\13\4\f\t\f\4\r\t\r\4\16\t\16\4\17\t\17\4\20\t\20\4\21\t\21\4\22"+
		"\t\22\4\23\t\23\4\24\t\24\3\2\3\2\3\3\3\3\3\4\3\4\3\5\3\5\3\6\3\6\3\7"+
		"\3\7\3\b\3\b\3\b\3\t\3\t\3\n\3\n\3\13\3\13\3\f\3\f\3\r\5\rB\n\r\3\r\3"+
		"\r\3\r\7\rG\n\r\f\r\16\rJ\13\r\3\r\3\r\5\rN\n\r\3\r\6\rQ\n\r\r\r\16\r"+
		"R\5\rU\n\r\3\16\5\16X\n\16\3\16\3\16\3\17\3\17\3\17\7\17_\n\17\f\17\16"+
		"\17b\13\17\5\17d\n\17\3\20\3\20\3\21\3\21\3\22\3\22\3\22\3\22\5\22n\n"+
		"\22\3\23\6\23q\n\23\r\23\16\23r\3\23\3\23\7\23w\n\23\f\23\16\23z\13\23"+
		"\3\23\3\23\3\23\7\23\177\n\23\f\23\16\23\u0082\13\23\3\23\5\23\u0085\n"+
		"\23\3\24\6\24\u0088\n\24\r\24\16\24\u0089\3\24\3\24\4x\u0080\2\25\3\3"+
		"\5\4\7\5\t\6\13\7\r\b\17\t\21\n\23\13\25\f\27\r\31\16\33\17\35\2\37\2"+
		"!\2#\20%\21\'\22\3\2\t\4\2GGgg\4\2--//\3\2\63;\3\2\62;\4\2JJTT\t\2\'("+
		",,\61;C\\aac|~~\5\2\13\f\17\17\"\"\u0098\2\3\3\2\2\2\2\5\3\2\2\2\2\7\3"+
		"\2\2\2\2\t\3\2\2\2\2\13\3\2\2\2\2\r\3\2\2\2\2\17\3\2\2\2\2\21\3\2\2\2"+
		"\2\23\3\2\2\2\2\25\3\2\2\2\2\27\3\2\2\2\2\31\3\2\2\2\2\33\3\2\2\2\2#\3"+
		"\2\2\2\2%\3\2\2\2\2\'\3\2\2\2\3)\3\2\2\2\5+\3\2\2\2\7-\3\2\2\2\t/\3\2"+
		"\2\2\13\61\3\2\2\2\r\63\3\2\2\2\17\65\3\2\2\2\218\3\2\2\2\23:\3\2\2\2"+
		"\25<\3\2\2\2\27>\3\2\2\2\31A\3\2\2\2\33W\3\2\2\2\35c\3\2\2\2\37e\3\2\2"+
		"\2!g\3\2\2\2#m\3\2\2\2%\u0084\3\2\2\2\'\u0087\3\2\2\2)*\7=\2\2*\4\3\2"+
		"\2\2+,\7*\2\2,\6\3\2\2\2-.\7.\2\2.\b\3\2\2\2/\60\7+\2\2\60\n\3\2\2\2\61"+
		"\62\7<\2\2\62\f\3\2\2\2\63\64\7%\2\2\64\16\3\2\2\2\65\66\7]\2\2\66\67"+
		"\7(\2\2\67\20\3\2\2\289\7_\2\29\22\3\2\2\2:;\7?\2\2;\24\3\2\2\2<=\7}\2"+
		"\2=\26\3\2\2\2>?\7\177\2\2?\30\3\2\2\2@B\7/\2\2A@\3\2\2\2AB\3\2\2\2BC"+
		"\3\2\2\2CD\5\35\17\2DH\7\60\2\2EG\5!\21\2FE\3\2\2\2GJ\3\2\2\2HF\3\2\2"+
		"\2HI\3\2\2\2IT\3\2\2\2JH\3\2\2\2KM\t\2\2\2LN\t\3\2\2ML\3\2\2\2MN\3\2\2"+
		"\2NP\3\2\2\2OQ\5!\21\2PO\3\2\2\2QR\3\2\2\2RP\3\2\2\2RS\3\2\2\2SU\3\2\2"+
		"\2TK\3\2\2\2TU\3\2\2\2U\32\3\2\2\2VX\7/\2\2WV\3\2\2\2WX\3\2\2\2XY\3\2"+
		"\2\2YZ\5\35\17\2Z\34\3\2\2\2[d\7\62\2\2\\`\5\37\20\2]_\5!\21\2^]\3\2\2"+
		"\2_b\3\2\2\2`^\3\2\2\2`a\3\2\2\2ad\3\2\2\2b`\3\2\2\2c[\3\2\2\2c\\\3\2"+
		"\2\2d\36\3\2\2\2ef\t\4\2\2f \3\2\2\2gh\t\5\2\2h\"\3\2\2\2in\t\6\2\2jk"+
		"\7N\2\2kl\7I\2\2ln\7V\2\2mi\3\2\2\2mj\3\2\2\2n$\3\2\2\2oq\t\7\2\2po\3"+
		"\2\2\2qr\3\2\2\2rp\3\2\2\2rs\3\2\2\2s\u0085\3\2\2\2tx\7$\2\2uw\13\2\2"+
		"\2vu\3\2\2\2wz\3\2\2\2xy\3\2\2\2xv\3\2\2\2y{\3\2\2\2zx\3\2\2\2{\u0085"+
		"\7$\2\2|\u0080\7)\2\2}\177\13\2\2\2~}\3\2\2\2\177\u0082\3\2\2\2\u0080"+
		"\u0081\3\2\2\2\u0080~\3\2\2\2\u0081\u0083\3\2\2\2\u0082\u0080\3\2\2\2"+
		"\u0083\u0085\7)\2\2\u0084p\3\2\2\2\u0084t\3\2\2\2\u0084|\3\2\2\2\u0085"+
		"&\3\2\2\2\u0086\u0088\t\b\2\2\u0087\u0086\3\2\2\2\u0088\u0089\3\2\2\2"+
		"\u0089\u0087\3\2\2\2\u0089\u008a\3\2\2\2\u008a\u008b\3\2\2\2\u008b\u008c"+
		"\b\24\2\2\u008c(\3\2\2\2\21\2AHMRTW`cmrx\u0080\u0084\u0089\3\b\2\2";
	public static final ATN _ATN =
		new ATNDeserializer().deserialize(_serializedATN.toCharArray());
	static {
		_decisionToDFA = new DFA[_ATN.getNumberOfDecisions()];
		for (int i = 0; i < _ATN.getNumberOfDecisions(); i++) {
			_decisionToDFA[i] = new DFA(_ATN.getDecisionState(i), i);
		}
	}
}