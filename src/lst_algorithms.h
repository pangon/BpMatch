/*

Copyright (C) 2003-2006 Christian Kreibich <christian@whoop.org>.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies of the Software and its documentation and acknowledgment shall be
given in the documentation and software packages that this Software was
used.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
#ifndef __lst_algorithms_h
#define __lst_algorithms_h

#include "lst_stree.h"

/**
 * LST_NodeVisitCB - callback signature for node visits.
 * @node: node currently visited.
 * @data: arbitrary data passed through.
 *
 * This is the signature of the callbacks used in several of
 * the algorithms below, that iterate over the tree. They
 * call a callback of this signature for every node they
 * visit.
 *
 * Returns: value > 0 if the iteration algorithm that called
 * this node is to proceed beyond this node, or 0 if not. Note
 * that this does not necessarily mean that an algorithm will
 * abort when the return value is 0.
 */
typedef int (*LST_NodeVisitCB)(LST_Node *node, void *data);


/**
 * lst_alg_bfs - breadth-first search of suffix tree.
 * @tree: suffix tree to iterate.
 * @callback: callback to call for each node.
 * @data: user data passed through to callback.
 *
 * The algorithm iterates the tree in breadth-first order, calling
 * @callback for each node visited.
 */
void            lst_alg_bfs(LST_STree *tree, LST_NodeVisitCB callback, void *data);


/**
 * lst_alg_dfs - depth-first search of suffix tree.
 * @tree: suffix tree to iterate.
 * @callback: callback to call for each node.
 * @data: user data passed through to callback.
 *
 * The algorithm iterates the tree in depth-first order, calling
 * @callback for each node visited.
 */
void            lst_alg_dfs(LST_STree *tree, LST_NodeVisitCB callback, void *data);


/**
 * lst_alg_bus - bottom-up search of suffix tree.
 * @tree: suffix tree to iterate.
 * @callback: callback to call for each node.
 * @data: user data passed through to callback.
 *
 * The algorithm iterates the tree in bottom-up order, calling @callback
 * for each node visited. This algorithm ignores the return value
 * of @callback.
 */
void            lst_alg_bus(LST_STree *tree, LST_NodeVisitCB callback, void *data);


/**
 * lst_alg_leafs - iterates all leafs in a suffix tree.
 * @tree: suffix tree to visit.
 * @callback: callback to call for each node.
 * @data: user data passed through to callback.
 *
 * The algorithm iterates over all leafs in the tree, calling @callback
 * for each node visited. If @callback returns 0, it stops.
 */
void            lst_alg_leafs(LST_STree *tree, LST_NodeVisitCB callback, void *data);


/**
 * lst_alg_set_visitors - builds visitor bitstrings in tree nodes.
 * @tree: tree to update.
 *
 * The algorithm updates the visitor elements in each node of @tree
 * to contain a one-bit for each string index that is contained
 * in the tree.
 *
 * Returns: bitstring representing a node visited by all strings.
 */
u_int           lst_alg_set_visitors(LST_STree *tree);


/**
 * lst_alg_longest_common_substring - computes the lcs for a suffix tree.
 * @tree: tree to use in computation.
 * @min_len: minimum length that common substrings must have to be returned.
 * @max_len: don't return strings longer than @max_len items.
 *
 * The algorithm computes the longest common substring(s) in @tree
 * and returns them as a new string set. This is currently a suboptimal
 * O(n^2) implementation until I have time for the more sophisticated
 * O(n) implementation available. If you want to limit the string length,
 * pass an appropriate value for @max_len, or pass 0 if you want the
 * longest string(s) possible. Similarly, if you want to receive only
 * longest common substrings of at least a certain number of items, use
 * @min_len for that, or pass 0 to indicate interest in everything.
 *
 * Returns: new string set, or %NULL when no strings were found.
 */
LST_StringSet  *lst_alg_longest_common_substring(LST_STree *tree,
						 u_int min_len,
						 u_int max_len);

/**
 * lst_alg_longest_repeated_substring - computes the lrs for a suffix tree.
 * @tree: tree to use in computation.
 * @min_len: minimum length that repeated substrings must have to be returned.
 * @max_len: don't return strings longer than @max_len items.
 *
 * The algorithm computes the longest repeated substring(s) in @tree
 * and returns them as a new string set. This is currently a suboptimal
 * O(n^2) implementation until I have time for the more sophisticated
 * O(n) implementation available. If you want to limit the string length,
 * pass an appropriate value for @max_len, or pass 0 if you want the
 * longest string(s) possible. Similarly, if you want to receive only
 * longest repeated substrings of at least a certain number of items, use
 * @min_len for that, or pass 0 to indicate interest in everything.
 *
 * Returns: new string set, or %NULL when no strings were found.
 */
LST_StringSet  *lst_alg_longest_repeated_substring(LST_STree *tree,
                                                  u_int min_len,
                                                  u_int max_len);

#endif
