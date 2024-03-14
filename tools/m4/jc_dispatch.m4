# Basic macro to detect support for libdispatch (aka Grand Central Dispatch)
# and compiler support for the block syntax extension commonly used with it

AC_DEFUN([JC_DISPATCH], [

AC_SEARCH_LIBS([dispatch_apply], [dispatch])
AC_CHECK_HEADER([dispatch/dispatch.h])

AC_MSG_CHECKING([whether a simple program can be compiled against libdispatch])
AC_LINK_IFELSE([AC_LANG_SOURCE([[

#include <stdlib.h>
#include <dispatch/dispatch.h>

static void kernel (void *context, size_t iteration)
{
    int *values = (int *) context;
    values[iteration] = iteration;
}

int main ()
{
    int *values = (int *) calloc(10, sizeof(int));
    dispatch_apply_f(10, DISPATCH_APPLY_AUTO, values, &kernel);
    free(values);
    return 0;
}

]])], [
    LIBDISPATCH_CPPFLAGS="-DHAVE_LIBDISPATCH"
    AC_MSG_RESULT([yes])
], [AC_MSG_RESULT([no])])


AC_MSG_CHECKING([compiler flag for block support])
orig_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[int main() { return ^{return 0;}(); }]])], [
    LIBDISPATCH_CPPFLAGS="$LIBDISPATCH_CPPFLAGS -DHAVE_BLOCKS"
    AC_MSG_RESULT([none required])
], [
    []_AC_LANG_PREFIX[]FLAGS="$orig_[]_AC_LANG_PREFIX[]FLAGS -fblocks"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([[int main() { return ^{return 0;}(); }]])], [
        LIBDISPATCH_CPPFLAGS="$LIBDISPATCH_CPPFLAGS -DHAVE_BLOCKS"
        LIBDISPATCH_[]_AC_LANG_PREFIX[]FLAGS="-fblocks"
        AC_MSG_RESULT([-fblocks])
    ], [
        AC_MSG_RESULT([none])
    ])
])
[]_AC_LANG_PREFIX[]FLAGS=$orig_[]_AC_LANG_PREFIX[]FLAGS

AC_SEARCH_LIBS([_Block_copy], [BlocksRuntime])
AC_CHECK_HEADER([Block.h])

])
