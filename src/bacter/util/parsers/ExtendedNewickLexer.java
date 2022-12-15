// Generated from /Users/vaughant/code/beast_and_friends/bacter/src/bacter/util/parsers/ExtendedNewick.g4 by ANTLR 4.10.1
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
	static { RuntimeMetaData.checkVersion("4.10.1", RuntimeMetaData.VERSION); }

	protected static final DFA[] _decisionToDFA;
	protected static final PredictionContextCache _sharedContextCache =
		new PredictionContextCache();
	public static final int
		T__0=1, T__1=2, T__2=3, T__3=4, T__4=5, T__5=6, T__6=7, T__7=8, T__8=9, 
		T__9=10, T__10=11, FLOAT=12, INT=13, RECTYPE=14, STRINGPRIM=15, WHITESPACE=16;
	public static String[] channelNames = {
		"DEFAULT_TOKEN_CHANNEL", "HIDDEN"
	};

	public static String[] modeNames = {
		"DEFAULT_MODE"
	};

	private static String[] makeRuleNames() {
		return new String[] {
			"T__0", "T__1", "T__2", "T__3", "T__4", "T__5", "T__6", "T__7", "T__8", 
			"T__9", "T__10", "FLOAT", "INT", "NNINT", "NZD", "D", "RECTYPE", "STRINGPRIM", 
			"WHITESPACE"
		};
	}
	public static final String[] ruleNames = makeRuleNames();

	private static String[] makeLiteralNames() {
		return new String[] {
			null, "';'", "'('", "','", "')'", "':'", "'#'", "'[&'", "']'", "'='", 
			"'{'", "'}'"
		};
	}
	private static final String[] _LITERAL_NAMES = makeLiteralNames();
	private static String[] makeSymbolicNames() {
		return new String[] {
			null, null, null, null, null, null, null, null, null, null, null, null, 
			"FLOAT", "INT", "RECTYPE", "STRINGPRIM", "WHITESPACE"
		};
	}
	private static final String[] _SYMBOLIC_NAMES = makeSymbolicNames();
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
	public String[] getChannelNames() { return channelNames; }

	@Override
	public String[] getModeNames() { return modeNames; }

	@Override
	public ATN getATN() { return _ATN; }

	public static final String _serializedATN =
		"\u0004\u0000\u0010\u008b\u0006\uffff\uffff\u0002\u0000\u0007\u0000\u0002"+
		"\u0001\u0007\u0001\u0002\u0002\u0007\u0002\u0002\u0003\u0007\u0003\u0002"+
		"\u0004\u0007\u0004\u0002\u0005\u0007\u0005\u0002\u0006\u0007\u0006\u0002"+
		"\u0007\u0007\u0007\u0002\b\u0007\b\u0002\t\u0007\t\u0002\n\u0007\n\u0002"+
		"\u000b\u0007\u000b\u0002\f\u0007\f\u0002\r\u0007\r\u0002\u000e\u0007\u000e"+
		"\u0002\u000f\u0007\u000f\u0002\u0010\u0007\u0010\u0002\u0011\u0007\u0011"+
		"\u0002\u0012\u0007\u0012\u0001\u0000\u0001\u0000\u0001\u0001\u0001\u0001"+
		"\u0001\u0002\u0001\u0002\u0001\u0003\u0001\u0003\u0001\u0004\u0001\u0004"+
		"\u0001\u0005\u0001\u0005\u0001\u0006\u0001\u0006\u0001\u0006\u0001\u0007"+
		"\u0001\u0007\u0001\b\u0001\b\u0001\t\u0001\t\u0001\n\u0001\n\u0001\u000b"+
		"\u0003\u000b@\b\u000b\u0001\u000b\u0001\u000b\u0001\u000b\u0005\u000b"+
		"E\b\u000b\n\u000b\f\u000bH\t\u000b\u0001\u000b\u0001\u000b\u0003\u000b"+
		"L\b\u000b\u0001\u000b\u0004\u000bO\b\u000b\u000b\u000b\f\u000bP\u0003"+
		"\u000bS\b\u000b\u0001\f\u0003\fV\b\f\u0001\f\u0001\f\u0001\r\u0001\r\u0001"+
		"\r\u0005\r]\b\r\n\r\f\r`\t\r\u0003\rb\b\r\u0001\u000e\u0001\u000e\u0001"+
		"\u000f\u0001\u000f\u0001\u0010\u0001\u0010\u0001\u0010\u0001\u0010\u0003"+
		"\u0010l\b\u0010\u0001\u0011\u0004\u0011o\b\u0011\u000b\u0011\f\u0011p"+
		"\u0001\u0011\u0001\u0011\u0005\u0011u\b\u0011\n\u0011\f\u0011x\t\u0011"+
		"\u0001\u0011\u0001\u0011\u0001\u0011\u0005\u0011}\b\u0011\n\u0011\f\u0011"+
		"\u0080\t\u0011\u0001\u0011\u0003\u0011\u0083\b\u0011\u0001\u0012\u0004"+
		"\u0012\u0086\b\u0012\u000b\u0012\f\u0012\u0087\u0001\u0012\u0001\u0012"+
		"\u0002v~\u0000\u0013\u0001\u0001\u0003\u0002\u0005\u0003\u0007\u0004\t"+
		"\u0005\u000b\u0006\r\u0007\u000f\b\u0011\t\u0013\n\u0015\u000b\u0017\f"+
		"\u0019\r\u001b\u0000\u001d\u0000\u001f\u0000!\u000e#\u000f%\u0010\u0001"+
		"\u0000\u0007\u0002\u0000EEee\u0002\u0000++--\u0001\u000019\u0001\u0000"+
		"09\u0002\u0000HHRR\u0007\u0000%&*+-9AZ__az||\u0003\u0000\t\n\r\r  \u0096"+
		"\u0000\u0001\u0001\u0000\u0000\u0000\u0000\u0003\u0001\u0000\u0000\u0000"+
		"\u0000\u0005\u0001\u0000\u0000\u0000\u0000\u0007\u0001\u0000\u0000\u0000"+
		"\u0000\t\u0001\u0000\u0000\u0000\u0000\u000b\u0001\u0000\u0000\u0000\u0000"+
		"\r\u0001\u0000\u0000\u0000\u0000\u000f\u0001\u0000\u0000\u0000\u0000\u0011"+
		"\u0001\u0000\u0000\u0000\u0000\u0013\u0001\u0000\u0000\u0000\u0000\u0015"+
		"\u0001\u0000\u0000\u0000\u0000\u0017\u0001\u0000\u0000\u0000\u0000\u0019"+
		"\u0001\u0000\u0000\u0000\u0000!\u0001\u0000\u0000\u0000\u0000#\u0001\u0000"+
		"\u0000\u0000\u0000%\u0001\u0000\u0000\u0000\u0001\'\u0001\u0000\u0000"+
		"\u0000\u0003)\u0001\u0000\u0000\u0000\u0005+\u0001\u0000\u0000\u0000\u0007"+
		"-\u0001\u0000\u0000\u0000\t/\u0001\u0000\u0000\u0000\u000b1\u0001\u0000"+
		"\u0000\u0000\r3\u0001\u0000\u0000\u0000\u000f6\u0001\u0000\u0000\u0000"+
		"\u00118\u0001\u0000\u0000\u0000\u0013:\u0001\u0000\u0000\u0000\u0015<"+
		"\u0001\u0000\u0000\u0000\u0017?\u0001\u0000\u0000\u0000\u0019U\u0001\u0000"+
		"\u0000\u0000\u001ba\u0001\u0000\u0000\u0000\u001dc\u0001\u0000\u0000\u0000"+
		"\u001fe\u0001\u0000\u0000\u0000!k\u0001\u0000\u0000\u0000#\u0082\u0001"+
		"\u0000\u0000\u0000%\u0085\u0001\u0000\u0000\u0000\'(\u0005;\u0000\u0000"+
		"(\u0002\u0001\u0000\u0000\u0000)*\u0005(\u0000\u0000*\u0004\u0001\u0000"+
		"\u0000\u0000+,\u0005,\u0000\u0000,\u0006\u0001\u0000\u0000\u0000-.\u0005"+
		")\u0000\u0000.\b\u0001\u0000\u0000\u0000/0\u0005:\u0000\u00000\n\u0001"+
		"\u0000\u0000\u000012\u0005#\u0000\u00002\f\u0001\u0000\u0000\u000034\u0005"+
		"[\u0000\u000045\u0005&\u0000\u00005\u000e\u0001\u0000\u0000\u000067\u0005"+
		"]\u0000\u00007\u0010\u0001\u0000\u0000\u000089\u0005=\u0000\u00009\u0012"+
		"\u0001\u0000\u0000\u0000:;\u0005{\u0000\u0000;\u0014\u0001\u0000\u0000"+
		"\u0000<=\u0005}\u0000\u0000=\u0016\u0001\u0000\u0000\u0000>@\u0005-\u0000"+
		"\u0000?>\u0001\u0000\u0000\u0000?@\u0001\u0000\u0000\u0000@A\u0001\u0000"+
		"\u0000\u0000AB\u0003\u001b\r\u0000BF\u0005.\u0000\u0000CE\u0003\u001f"+
		"\u000f\u0000DC\u0001\u0000\u0000\u0000EH\u0001\u0000\u0000\u0000FD\u0001"+
		"\u0000\u0000\u0000FG\u0001\u0000\u0000\u0000GR\u0001\u0000\u0000\u0000"+
		"HF\u0001\u0000\u0000\u0000IK\u0007\u0000\u0000\u0000JL\u0007\u0001\u0000"+
		"\u0000KJ\u0001\u0000\u0000\u0000KL\u0001\u0000\u0000\u0000LN\u0001\u0000"+
		"\u0000\u0000MO\u0003\u001f\u000f\u0000NM\u0001\u0000\u0000\u0000OP\u0001"+
		"\u0000\u0000\u0000PN\u0001\u0000\u0000\u0000PQ\u0001\u0000\u0000\u0000"+
		"QS\u0001\u0000\u0000\u0000RI\u0001\u0000\u0000\u0000RS\u0001\u0000\u0000"+
		"\u0000S\u0018\u0001\u0000\u0000\u0000TV\u0005-\u0000\u0000UT\u0001\u0000"+
		"\u0000\u0000UV\u0001\u0000\u0000\u0000VW\u0001\u0000\u0000\u0000WX\u0003"+
		"\u001b\r\u0000X\u001a\u0001\u0000\u0000\u0000Yb\u00050\u0000\u0000Z^\u0003"+
		"\u001d\u000e\u0000[]\u0003\u001f\u000f\u0000\\[\u0001\u0000\u0000\u0000"+
		"]`\u0001\u0000\u0000\u0000^\\\u0001\u0000\u0000\u0000^_\u0001\u0000\u0000"+
		"\u0000_b\u0001\u0000\u0000\u0000`^\u0001\u0000\u0000\u0000aY\u0001\u0000"+
		"\u0000\u0000aZ\u0001\u0000\u0000\u0000b\u001c\u0001\u0000\u0000\u0000"+
		"cd\u0007\u0002\u0000\u0000d\u001e\u0001\u0000\u0000\u0000ef\u0007\u0003"+
		"\u0000\u0000f \u0001\u0000\u0000\u0000gl\u0007\u0004\u0000\u0000hi\u0005"+
		"L\u0000\u0000ij\u0005G\u0000\u0000jl\u0005T\u0000\u0000kg\u0001\u0000"+
		"\u0000\u0000kh\u0001\u0000\u0000\u0000l\"\u0001\u0000\u0000\u0000mo\u0007"+
		"\u0005\u0000\u0000nm\u0001\u0000\u0000\u0000op\u0001\u0000\u0000\u0000"+
		"pn\u0001\u0000\u0000\u0000pq\u0001\u0000\u0000\u0000q\u0083\u0001\u0000"+
		"\u0000\u0000rv\u0005\"\u0000\u0000su\t\u0000\u0000\u0000ts\u0001\u0000"+
		"\u0000\u0000ux\u0001\u0000\u0000\u0000vw\u0001\u0000\u0000\u0000vt\u0001"+
		"\u0000\u0000\u0000wy\u0001\u0000\u0000\u0000xv\u0001\u0000\u0000\u0000"+
		"y\u0083\u0005\"\u0000\u0000z~\u0005\'\u0000\u0000{}\t\u0000\u0000\u0000"+
		"|{\u0001\u0000\u0000\u0000}\u0080\u0001\u0000\u0000\u0000~\u007f\u0001"+
		"\u0000\u0000\u0000~|\u0001\u0000\u0000\u0000\u007f\u0081\u0001\u0000\u0000"+
		"\u0000\u0080~\u0001\u0000\u0000\u0000\u0081\u0083\u0005\'\u0000\u0000"+
		"\u0082n\u0001\u0000\u0000\u0000\u0082r\u0001\u0000\u0000\u0000\u0082z"+
		"\u0001\u0000\u0000\u0000\u0083$\u0001\u0000\u0000\u0000\u0084\u0086\u0007"+
		"\u0006\u0000\u0000\u0085\u0084\u0001\u0000\u0000\u0000\u0086\u0087\u0001"+
		"\u0000\u0000\u0000\u0087\u0085\u0001\u0000\u0000\u0000\u0087\u0088\u0001"+
		"\u0000\u0000\u0000\u0088\u0089\u0001\u0000\u0000\u0000\u0089\u008a\u0006"+
		"\u0012\u0000\u0000\u008a&\u0001\u0000\u0000\u0000\u000f\u0000?FKPRU^a"+
		"kpv~\u0082\u0087\u0001\u0006\u0000\u0000";
	public static final ATN _ATN =
		new ATNDeserializer().deserialize(_serializedATN.toCharArray());
	static {
		_decisionToDFA = new DFA[_ATN.getNumberOfDecisions()];
		for (int i = 0; i < _ATN.getNumberOfDecisions(); i++) {
			_decisionToDFA[i] = new DFA(_ATN.getDecisionState(i), i);
		}
	}
}