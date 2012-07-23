/*

Copyright (C) 2000 - 2006 Christian Kreibich <christian@whoop.org>.

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
#ifndef __lst_debug_h
#define __lst_debug_h

#include "lst_stree.h"


/**
 * lst_debug_print_tree - prints out a tree to stderr.
 * @tree: tree to print.
 *
 * We make this always available, not only with -DDEBUG.
 */
void   lst_debug_print_tree(LST_STree *tree);



#ifdef LST_DEBUG

/**
 * D_PRINT_TREE - prints out a tree to stderr.
 * @tree: tree to print.
 *
 * This is just a wrapper around lst_debug_print_tree(),
 * so nothing will be printed if debugging isn't enabled.
 */
#define D_PRINT_TREE(tree) lst_debug_print_tree(tree)

char  *lst_debug_print_substring(LST_String *string,
				 u_int start_index,
				 u_int end_index,
				 u_int extra_index);

#define lst_stderr(...) fprintf(stderr, __VA_ARGS__)

/**
 * D - prints debugging output
 * @x: debugging information.
 *
 * Use this macro to output debugging information. @x is
 * the content as you would pass it to printf(), including
 * braces to make the arguments appear as one argument to
 * the macro. The macro is automatically deleted if -DDEBUG
 * is not passed at build time.
 */
#define D(x)                  do { fprintf(stderr, "%s/%i: ", __FILE__, __LINE__); lst_stderr x ; } while (0)

/**
 * D_ASSERT - debugging assertion.
 * @exp: expression to evaluate.
 * @msg: message to output if @exp fails.
 *
 * The macro outputs @msg if the expression @exp evaluates
 * to %FALSE.
 */
#define D_ASSERT(exp, msg)    if (! exp) { fprintf(stderr, "%s/%i: %s\n", __FILE__, __LINE__, msg); }

/**
 * D_ASSERT_PTR - pointer existence assertion.
 * @ptr: pointer to check.
 *
 * The macro asserts the existence (i.e. non-NULL-ness) of
 * the given pointer, and outpus a message if it is %NULL.
 */
#define D_ASSERT_PTR(ptr)     D_ASSERT(ptr, "pointer is NULL.")

#else
#define D_PRINT_TREE(x)
#define lst_debug_print_substring(a, b, c, d) NULL
#define D(x)                  
#define D_ASSERT(exp, msg)    
#define D_ASSERT_PTR(ptr)     
#endif

#endif 

