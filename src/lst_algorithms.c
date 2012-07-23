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
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lst_string.h"
#include "lst_stree.h"
#include "lst_structs.h"
#include "lst_debug.h"
#include "lst_algorithms.h"

typedef struct lst_node_it
{
  TAILQ_ENTRY(lst_node_it) items;
  LST_Node *node;
} LST_NodeIt;


static LST_NodeIt *
alg_node_it_new(LST_Node *node)
{
  LST_NodeIt *it;

  it = (LST_NodeIt*)calloc(1, sizeof(LST_NodeIt));
  it->node = node;

  return it;
}

static void
alg_node_it_free(LST_NodeIt *it)
{
  if (!it)
    return;

  free(it);
}

void
lst_alg_bfs(LST_STree *tree, LST_NodeVisitCB callback, void *data)
{
  LST_Node *node;
  LST_Edge *edge;
  TAILQ_HEAD(qhead, lst_node) queue;

  if (!tree || !callback)
    return;

  TAILQ_INIT(&queue);
  TAILQ_INSERT_HEAD(&queue, tree->root_node, iteration);

  while (queue.tqh_first)
    {
      node = queue.tqh_first;
      TAILQ_REMOVE(&queue, queue.tqh_first, iteration);

      if (callback(node, data))
	{
	  for (edge = node->kids.lh_first; edge; edge = edge->siblings.le_next)
	    TAILQ_INSERT_TAIL(&queue, edge->dst_node, iteration);
	}
    }

}


void
lst_alg_dfs(LST_STree *tree, LST_NodeVisitCB callback, void *data)
{
  LST_Node *node;
  LST_Edge *edge;
  TAILQ_HEAD(shead, lst_node) stack;

  if (!tree || !callback)
    return;

  TAILQ_INIT(&stack);
  TAILQ_INSERT_HEAD(&stack, tree->root_node, iteration);

  while (stack.tqh_first)
    {
      node = stack.tqh_first;
      TAILQ_REMOVE(&stack, stack.tqh_first, iteration);

      if (callback(node, data))
	{
	  for (edge = node->kids.lh_first; edge; edge = edge->siblings.le_next)
	    TAILQ_INSERT_HEAD(&stack, edge->dst_node, iteration);
	}
    }
}


static int
alg_clear_busflag(LST_Node *node, void *data)
{
  node->bus_visited = 0;

  return 1;
  data = NULL;
}


void
lst_alg_bus(LST_STree *tree, LST_NodeVisitCB callback, void *data)
{
  TAILQ_HEAD(nodes_s, lst_node_it) nodes;
  LST_Node *node;
  LST_NodeIt *it;

  TAILQ_INIT(&nodes);
  lst_alg_bfs(tree, alg_clear_busflag, NULL);

  if (!tree->leafs.lh_first)
    return;

  for (node = tree->leafs.lh_first; node; node = node->leafs.le_next)
    {
      callback(node, data);

      if (node == tree->root_node)
	continue;

      node->up_edge->src_node->bus_visited++;

      if (node->up_edge->src_node->bus_visited == 1)
	{
	  it = alg_node_it_new(node->up_edge->src_node);
	  TAILQ_INSERT_TAIL(&nodes, it, items);
	}
    }

  while (nodes.tqh_first)
    {
      it = nodes.tqh_first;
      node = it->node;
      TAILQ_REMOVE(&nodes, nodes.tqh_first, items);

      if (node->bus_visited < node->num_kids)
	{
	  TAILQ_INSERT_TAIL(&nodes, it, items);
	  continue;
	}

      callback(node, data);
      alg_node_it_free(it);

      if (node == tree->root_node)
        continue;

      node->up_edge->src_node->bus_visited++;

      if (node->up_edge->src_node->bus_visited == 1)
	{
	  it = alg_node_it_new(node->up_edge->src_node);
	  TAILQ_INSERT_TAIL(&nodes, it, items);
	}
    }
}


void
lst_alg_leafs(LST_STree *tree, LST_NodeVisitCB callback, void *data)
{
  LST_Node *node;

  if (!tree || !callback)
    return;

  for (node = tree->leafs.lh_first; node; node = node->leafs.le_next)
    {
      if (callback(node, data) == 0)
	break;
    }
}


typedef struct lst_lcs_data
{
  LST_STree    *tree;

  /* 1 for longest common substring, 0 for longest repeated substring */
  int          lcs;

  /* Used in LCS case only */
  u_int         all_visitors;

  TAILQ_HEAD(nodes, lst_node_it) nodes;
  int           deepest;
  int           num_deepest;
  int           max_depth;

} LST_LCS_Data;

static int
alg_clear_visitors(LST_Node *node, void *data)
{
  node->visitors = 0;

  return 1;
  data = NULL;
}

static int
alg_set_visitors(LST_Node *node, LST_LCS_Data *data)
{
  if (lst_node_is_root(node))
    {
      D(("Node %u: visitors %i\n", node->id, node->visitors));
      return 1;
    }

  if (lst_node_is_leaf(node))
    {
      int index =
	1 << lst_stree_get_string_index(data->tree, node->up_edge->range.string);

      node->visitors = index;
      node->up_edge->src_node->visitors |= index;
    }
  else
    {
      node->up_edge->src_node->visitors |= node->visitors;
    }

  if (node->up_edge->src_node->visitors > data->all_visitors)
    data->all_visitors = node->up_edge->src_node->visitors;

  D(("Node %u: visitors %i\n", node->id, node->visitors));
  return 1;
}

u_int
lst_alg_set_visitors(LST_STree *tree)
{
  LST_LCS_Data data;

  if (!tree)
    return 0;

  if (!tree->needs_visitor_update)
    return tree->visitors;

  memset(&data, 0, sizeof(LST_LCS_Data));
  data.tree = tree;

  /* First, establish the visitor bitstrings in the tree. */
  lst_alg_bus(tree, alg_clear_visitors, NULL);
  lst_alg_bus(tree, (LST_NodeVisitCB) alg_set_visitors, &data);

  tree->needs_visitor_update = 0;
  tree->visitors = data.all_visitors;

  return data.all_visitors;
}


static int
alg_find_deepest(LST_Node *node, LST_LCS_Data *data)
{
  LST_NodeIt *it;
  int depth = lst_node_get_string_length(node);

  if (data->lcs)
    {
      if (node->visitors != data->all_visitors)
       return 0;
    }
  else
    {
      if (node->num_kids < 1)
       return 0;
    }

  if (data->deepest <= data->max_depth)
    {
      if (depth >= data->deepest)
	{
	  it = alg_node_it_new(node);

	  if (depth > data->deepest)
	    {
	      data->deepest = depth;
	      data->num_deepest = 0;
	    }

	  data->num_deepest++;
	  TAILQ_INSERT_HEAD(&data->nodes, it, items);
	}
    }
  else if (depth >= data->max_depth)
    {
      it = alg_node_it_new(node);
      data->num_deepest++;
      TAILQ_INSERT_HEAD(&data->nodes, it, items);
    }

  return 1;
}

static LST_StringSet *
alg_longest_substring(LST_STree *tree, u_int min_len, u_int max_len, int lcs)
{
  LST_StringSet *result = NULL;
  LST_String *string;
  LST_LCS_Data data;
  LST_NodeIt *it;

  if (!tree)
    return NULL;

  memset(&data, 0, sizeof(LST_LCS_Data));
  data.tree = tree;
  data.lcs = lcs;
  if (lcs)
    data.all_visitors = lst_alg_set_visitors(tree);

  if (max_len > 0)
    data.max_depth = (int) max_len;
  else
    data.max_depth = INT_MAX;

  TAILQ_INIT(&data.nodes);

  /* Now do a DSF finding the node with the largest string-
   * depth that has all strings as visitors.
   */
  lst_alg_dfs(tree, (LST_NodeVisitCB) alg_find_deepest, &data);
  D(("Deepest nodes found -- we have %u longest substring(s) at depth %u.\n",
     data.num_deepest, data.deepest));

  /* Now, data.num_deepest tells us how many largest substrings
   * we have, and the first num_deepest items in data.nodes are
   * the end nodes in the suffix tree that define these substrings.
   */
  while ( (it = data.nodes.tqh_first))
    {
      if (--data.num_deepest >= 0)
	{
	  /* Get our longest common string's length, and if it's
	   * long enough for our requirements, put it in the result
	   * set. We need to allocate that first if we haven't yet
	   * inserted any strings.
	   */
	  if ((u_int) lst_node_get_string_length(it->node) >= min_len)
	    {
	      string = lst_node_get_string(it->node, (int) max_len);

	      if (!result)
		result = lst_stringset_new();

	      lst_stringset_add(result, string);
	    }
	}

      TAILQ_REMOVE(&data.nodes, it, items);
      alg_node_it_free(it);
    }

  return result;
}

LST_StringSet *
lst_alg_longest_common_substring(LST_STree *tree, u_int min_len, u_int max_len)
{
  return alg_longest_substring(tree, min_len, max_len, 1);
}

LST_StringSet *
lst_alg_longest_repeated_substring(LST_STree *tree, u_int min_len, u_int max_len)
{
  return alg_longest_substring(tree, min_len, max_len, 0);
}

