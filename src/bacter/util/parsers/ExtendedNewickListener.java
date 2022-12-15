// Generated from /Users/vaughant/code/beast_and_friends/bacter/src/bacter/util/parsers/ExtendedNewick.g4 by ANTLR 4.10.1
package bacter.util.parsers;
import org.antlr.v4.runtime.tree.ParseTreeListener;

/**
 * This interface defines a complete listener for a parse tree produced by
 * {@link ExtendedNewickParser}.
 */
public interface ExtendedNewickListener extends ParseTreeListener {
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#tree}.
	 * @param ctx the parse tree
	 */
	void enterTree(ExtendedNewickParser.TreeContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#tree}.
	 * @param ctx the parse tree
	 */
	void exitTree(ExtendedNewickParser.TreeContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#node}.
	 * @param ctx the parse tree
	 */
	void enterNode(ExtendedNewickParser.NodeContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#node}.
	 * @param ctx the parse tree
	 */
	void exitNode(ExtendedNewickParser.NodeContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#post}.
	 * @param ctx the parse tree
	 */
	void enterPost(ExtendedNewickParser.PostContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#post}.
	 * @param ctx the parse tree
	 */
	void exitPost(ExtendedNewickParser.PostContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#label}.
	 * @param ctx the parse tree
	 */
	void enterLabel(ExtendedNewickParser.LabelContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#label}.
	 * @param ctx the parse tree
	 */
	void exitLabel(ExtendedNewickParser.LabelContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#hybrid}.
	 * @param ctx the parse tree
	 */
	void enterHybrid(ExtendedNewickParser.HybridContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#hybrid}.
	 * @param ctx the parse tree
	 */
	void exitHybrid(ExtendedNewickParser.HybridContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#meta}.
	 * @param ctx the parse tree
	 */
	void enterMeta(ExtendedNewickParser.MetaContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#meta}.
	 * @param ctx the parse tree
	 */
	void exitMeta(ExtendedNewickParser.MetaContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#attrib}.
	 * @param ctx the parse tree
	 */
	void enterAttrib(ExtendedNewickParser.AttribContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#attrib}.
	 * @param ctx the parse tree
	 */
	void exitAttrib(ExtendedNewickParser.AttribContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#attribValue}.
	 * @param ctx the parse tree
	 */
	void enterAttribValue(ExtendedNewickParser.AttribValueContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#attribValue}.
	 * @param ctx the parse tree
	 */
	void exitAttribValue(ExtendedNewickParser.AttribValueContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#number}.
	 * @param ctx the parse tree
	 */
	void enterNumber(ExtendedNewickParser.NumberContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#number}.
	 * @param ctx the parse tree
	 */
	void exitNumber(ExtendedNewickParser.NumberContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#vector}.
	 * @param ctx the parse tree
	 */
	void enterVector(ExtendedNewickParser.VectorContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#vector}.
	 * @param ctx the parse tree
	 */
	void exitVector(ExtendedNewickParser.VectorContext ctx);
	/**
	 * Enter a parse tree produced by {@link ExtendedNewickParser#string}.
	 * @param ctx the parse tree
	 */
	void enterString(ExtendedNewickParser.StringContext ctx);
	/**
	 * Exit a parse tree produced by {@link ExtendedNewickParser#string}.
	 * @param ctx the parse tree
	 */
	void exitString(ExtendedNewickParser.StringContext ctx);
}