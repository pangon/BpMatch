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

#include <sys/param.h>

#include "lst_string.h"
#include "lst_macros.h"
#include "lst_debug.h"

static int   string_id_counter;

static int   string_byte_cmp_func(char *item1, char *item2);
static void  string_byte_copy_func(char *src, char *dst);
static char *string_print_func(LST_StringIndex *index);

static LST_StringClass byte_class =
  {
    (LST_StringItemCmpFunc) string_byte_cmp_func,
    (LST_StringItemCopyFunc) string_byte_copy_func,
    string_print_func
  };

static int
string_byte_cmp_func(char *item1, char *item2)
{
  if (*item1 < *item2)
    return -1;

  if (*item1 > *item2)
    return 1;

  return 0;
}


static void
string_byte_copy_func(char *src, char *dst)
{
  *dst = *src;
}


static char *
string_print_func(LST_StringIndex *index)
{
  static char s[3][4096];
  static int  i = 0;
  int idx;

  if (index->start_index == index->string->num_items - 1)
    return "<eos>";

  memcpy(s[i], ((char *) index->string->data) + index->start_index,
	 *(index->end_index) - index->start_index + 1);
  s[i][*(index->end_index) - index->start_index + 1] = '\0';

  if (index->extra_index > 0)
    {
      char tmp[128];

      snprintf(tmp, 128, "[%c]", *(((char *) index->string->data) + index->extra_index));
      strcat(s[i], tmp);
    }

  idx = i;
  i = (i + 1) % 3;

  return s[idx];
}


LST_String    *
lst_string_new(void *data, u_int item_size, u_int num_items)
{
  LST_String *string;

  if (item_size == 0 || num_items == 0)
    return NULL;

  string = (LST_String*)calloc(1, sizeof(LST_String));
  if (!string)
    {
      D(("Out of memory.\n"));
      return NULL;
    }

  string->id = ++string_id_counter;

  /* Logically, we want one more item than given; we treat that as a
   * special end-of-string marker so that no suffix of our string can ever
   * be the prefix of another suffix. For the problems that this
   * would cause, see Gusfield.
   */
  string->num_items  = num_items + 1;
  string->item_size  = item_size;
  string->sclass     = &byte_class;

  string->data = calloc(num_items, item_size);

  if (!string->data)
    {
      D(("Out of memory.\n"));
      free(string);
      return NULL;
    }

  if (data)
    memcpy(string->data, data, item_size * num_items);

  return string;
}


void
lst_string_init(LST_String *string, void *data, u_int item_size, u_int num_items)
{
  if (!string || !data || item_size == 0)
    {
      D(("String not initialized!\n"));
      return;
    }

  memset(string, 0, sizeof(LST_String));

  string->id = ++string_id_counter;

  string->data          = data;
  string->data_external = 1;
  string->num_items     = num_items + 1;
  string->item_size     = item_size;
  string->sclass        = &byte_class;
}


void *
lst_string_get_item(LST_String *string, u_int index)
{
  char *data = (char *) string->data;

  if (index >= string->num_items)
    return NULL;

  return (void *)(data + index * string->item_size);
}


void
lst_string_free(LST_String *string)
{
  if (!string)
    return;

  if (string->data && !string->data_external)
    free(string->data);

  free(string);
}


void *
lst_string_free_keep_data(LST_String *string)
{
  void *data;

  if (!string)
    return NULL;

  data = string->data;
  free(string);

  return data;
}


u_int
lst_string_get_length(LST_String *string)
{
  if (!string)
    return 0;

  return string->num_items - 1;
}


const char    *
lst_string_print(LST_String *string)
{
  LST_StringIndex tmp_range;

  if (!string)
    return NULL;

  lst_string_index_init(&tmp_range);

  tmp_range.string = string;
  tmp_range.start_index  = 0;
  *(tmp_range.end_index) = string->num_items - 1;
  tmp_range.extra_index  = 0;

  return string->sclass->print_func(&tmp_range);
}


void
lst_string_item_copy(LST_String *src, u_int src_index,
		     LST_String *dst, u_int dst_index)
{
  void *src_item, *dst_item;

  if (!src || !dst || src_index >= src->num_items || dst_index >= dst->num_items)
    return;

  src_item = lst_string_get_item(src, src_index);
  dst_item = lst_string_get_item(dst, dst_index);

  src->sclass->copy_func(src_item, dst_item);
}


int
lst_string_eq(LST_String *s1, u_int item1,
	      LST_String *s2, u_int item2)
{
  if (!s1 || !s2 || item1 >= s1->num_items || item2 >= s2->num_items)
    return 0;

  /* Treat the end-of-string markers separately: */
  if (item1 == s1->num_items - 1 || item2 == s2->num_items - 1) {
    if (item1 == s1->num_items - 1 && item2 == s2->num_items - 1) {
      if (s1 == s2)
	{
	  D(("Comparing end of identical strings\n"));
	  return 1;
	} else {
	  D(("Comparing end of different strings\n"));
	  return 0;
        }
     } else {
    	D(("Comparing end and non-end\n"));
        return 0;
     }
  }

  return !(s1->sclass->cmp_func(lst_string_get_item(s1, item1),
				lst_string_get_item(s2, item2)));
}


u_int
lst_string_items_common(LST_String *s1, u_int off1,
			LST_String *s2, u_int off2,
			u_int max_len)
{
  u_int i = 0;
  u_int len;

  if (!s1 || !s2 || off1 >= s1->num_items || off2 >= s2->num_items)
    return 0;

  len = MIN(MIN(s1->num_items - off1, s2->num_items - off2), max_len);

  for (i = 0; i < len; i++)
    {
      if (! lst_string_eq(s1, off1 + i,
			  s2, off2 + i))
	return i;
    }

  return len;
}


void
lst_string_index_init(LST_StringIndex *index)
{
  if (!index)
    return;

  memset(index, 0, sizeof(LST_StringIndex));
  index->end_index = &index->end_index_local;
}


void
lst_string_index_copy(LST_StringIndex *src, LST_StringIndex *dst)
{
  if (!src || !dst)
    return;

  dst->string = src->string;
  dst->start_index = src->start_index;
  *(dst->end_index) = *(src->end_index);
  dst->extra_index = src->extra_index;
}


char *
lst_string_print_hex(LST_StringIndex *index)
{
  char *s;
  int tmp;
  u_int i, j, end_index, len;
  u_char *data, *data_end;
  char *sptr;

  end_index = *(index->end_index);
  if (end_index == index->string->num_items - 1)
    end_index--;

  len = end_index - index->start_index + 1;

  if (index->start_index == index->string->num_items - 1)
    return "<eos>";

  /* Allocate a string the can hold the string data plus whitespace etc */
  s = (char*)calloc((2 * len) + 18 * len / 16 + 10, sizeof(char));
  if (!s)
    return NULL;

  sptr = s;

  if (index->start_index != LST_EMPTY_STRING)
    {
      data = (u_char*)(index->string->data) + index->start_index;
      data_end = data + len;

      for (i = 0; i < len; i += 16)
	{
	  for (j = 0; j < 16 && data < data_end; j++, data++)
	    {
	      sprintf(sptr, "%.2X ", *data);
	      sptr += 3;
	    }

	  /* Insert a newline if we have more to print */
	  if (i + 16 < len)
	    {
	      *sptr = '\n';
	      sptr++;
	    }
	}
    }

  if (index->extra_index > 0)
    {
      sprintf(sptr, "[%.2X]", *(((u_char *) index->string->data) + index->extra_index));
      sptr += 4;
    }

  *sptr = '\0';

  return s;
}


LST_StringClass *
lst_stringclass_new(LST_StringItemCmpFunc cmp_func,
		    LST_StringItemCopyFunc copy_func,
		    LST_StringPrintFunc print_func)
{
  LST_StringClass *sclass;

  sclass = (LST_StringClass*)calloc(1, sizeof(LST_StringClass));
  if (!sclass)
    return NULL;

  sclass->cmp_func   = cmp_func ? cmp_func : (LST_StringItemCmpFunc) string_byte_cmp_func;
  sclass->copy_func  = copy_func ? copy_func : (LST_StringItemCopyFunc) string_byte_copy_func;
  sclass->print_func = print_func ? print_func : string_print_func;

  return sclass;
}


void
lst_stringclass_free(LST_StringClass *sclass)
{
  if (sclass)
    free(sclass);
}


void
lst_stringclass_set_defaults(LST_StringItemCmpFunc cmp_func,
			     LST_StringItemCopyFunc copy_func,
			     LST_StringPrintFunc print_func)
{
  byte_class.cmp_func   = cmp_func ? cmp_func : (LST_StringItemCmpFunc) string_byte_cmp_func;
  byte_class.copy_func  = copy_func ? copy_func : (LST_StringItemCopyFunc) string_byte_copy_func;
  byte_class.print_func = print_func ? print_func : string_print_func;
}


LST_StringClass *
lst_string_set_class(LST_String *string, LST_StringClass *sclass)
{
  LST_StringClass *old;

  if (!string)
    return NULL;

  old = string->sclass;
  string->sclass = sclass ? sclass : &byte_class;

  return old;
}


LST_StringSet *
lst_stringset_new(void)
{
  LST_StringSet *set;

  set = (LST_StringSet*)calloc(1, sizeof(LST_StringSet));

  if (!set)
    return NULL;

  LIST_INIT(&set->members);

  return set;
}


void
lst_stringset_add(LST_StringSet *set, LST_String *string)
{
  if (!set || !string)
    return;

  LIST_INSERT_HEAD(&set->members, string, set);
  set->size++;
}


void
lst_stringset_remove(LST_StringSet *set, LST_String *string)
{
  LST_String *set_string;

  if (!set || !string)
    return;

  for (set_string = set->members.lh_first; set_string; set_string = set_string->set.le_next)
    {
      if (set_string->id != string->id)
	continue;

      LIST_REMOVE(string, set);
      set->size--;
      return;
    }
}


void
lst_stringset_foreach(LST_StringSet *set,
		      LST_StringCB callback,
		      void *user_data)
{
  LST_String *string;

  if (!set || !callback)
    return;

  for (string = set->members.lh_first; string; string = string->set.le_next)
    callback(string, user_data);
}


void
lst_stringset_free(LST_StringSet *set)
{
  LST_String *string;

  if (!set)
    return;

  while (set->members.lh_first)
    {
      string = set->members.lh_first;
      LIST_REMOVE(set->members.lh_first, set);
      lst_string_free(string);
    }

  free(set);
}
