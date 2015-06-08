// Generated from /home/tvaughan/code/beast_and_friends/bacter/src/bacter/util/parsers/ExtendedNewick.g4 by ANTLR 4.5
package bacter.util.parsers;
import org.antlr.v4.runtime.misc.NotNull;
import org.antlr.v4.runtime.tree.ParseTreeVisitor;

/**
 * This interface defines a complete generic visitor for a parse tree produced
 * by {@link ExtendedNewickParser}.
 *
 * @param <T> The return type of the visit operation. Use {@link Void} for
 * operations with no return type.
 */
public interface ExtendedNewickVisitor<T> extends ParseTreeVisitor<T> {
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#tree}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitTree(@NotNull ExtendedNewickParser.TreeContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#node}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitNode(@NotNull ExtendedNewickParser.NodeContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#post}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitPost(@NotNull ExtendedNewickParser.PostContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#label}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitLabel(@NotNull ExtendedNewickParser.LabelContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#hybrid}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitHybrid(@NotNull ExtendedNewickParser.HybridContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#meta}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitMeta(@NotNull ExtendedNewickParser.MetaContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#attrib}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitAttrib(@NotNull ExtendedNewickParser.AttribContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#attribValue}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitAttribValue(@NotNull ExtendedNewickParser.AttribValueContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#number}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitNumber(@NotNull ExtendedNewickParser.NumberContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#vector}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitVector(@NotNull ExtendedNewickParser.VectorContext ctx);
	/**
	 * Visit a parse tree produced by {@link ExtendedNewickParser#string}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitString(@NotNull ExtendedNewickParser.StringContext ctx);
}