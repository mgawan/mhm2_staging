/*
 * Tessellation Many-core Operating System. V2.
 *
 * Cell runtime.
 *
 * Copyright (c) 2014 The Regents of the University of California.
 *
 * This code is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 2 only, as published by
 * the Free Software Foundation.
 *
 * This code is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 2 only, as published by
 * the Free Software Foundation.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU General Public License version 2 for more
 * details (a copy is included in the LICENSE file that accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version 2
 * along with this work; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#ifndef __COLORS_H
#define __COLORS_H

#include <array>

#ifdef USE_COLORS

#define KNORM  "\x1B[0m"
#define KBLACK "\x1B[30m"
#define KRED  "\x1B[31m"
#define KGREEN  "\x1B[32m"
#define KYELLOW  "\x1B[33m"
#define KBLUE  "\x1B[34m"
#define KMAGENTA  "\x1B[35m"
#define KCYAN  "\x1B[36m"
#define KWHITE "\x1B[37m"
#define KGRAY "\x1B[90m"
#define KLRED "\x1B[91m"
#define KLGREEN "\x1B[92m"
#define KLYELLOW "\x1B[93m"
#define KLBLUE "\x1B[94m"
#define KLMAGENTA "\x1B[95m"
#define KLCYAN "\x1B[96m"
#define KLWHITE  "\x1B[97m"
#define KWHITEBK "\x1B[107m"
#define KGREENBK "\x1B[44m"
#define KLGREENBK "\x1B[102m"

static const std::array<std::string, 17> COLORS = {KNORM, KBLACK, KRED, KGREEN, KYELLOW, KBLUE, KMAGENTA, KCYAN, KWHITE, KGRAY,
                                                   KLRED, KLGREEN, KLYELLOW, KLBLUE, KLMAGENTA, KLCYAN, KLWHITE};
                                                  
#else

#define KNORM  ""
#define KBLACK ""
#define KRED  ""
#define KGREEN  ""
#define KYELLOW  ""
#define KBLUE  ""
#define KMAGENTA  ""
#define KCYAN  ""
#define KWHITE ""
#define KGRAY ""
#define KLRED ""
#define KLGREEN ""
#define KLYELLOW ""
#define KLBLUE ""
#define KLMAGENTA ""
#define KLCYAN ""
#define KLWHITE  ""
#define KWHITEBK ""
#define KGREENBK ""
#define KLGREENBK ""
#endif

#endif
