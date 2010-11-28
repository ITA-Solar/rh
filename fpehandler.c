/* ------- file: -------------------------- fpehandler.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Oct 15 16:26:36 2002 --

       --------------------------                      ----------RH-- */

/* --- Trap floating point exceptions on various machines:

         SGI mips, IRIX 5.2      (compile with: -Dmips -DIRIX5)
         SUN sparc, SunOS 4.x    (compile with: -Dsparc -DSunOS4)
         SUN sparc, SunOS 5.x    (compile with: -Dsparc -DSunOS5)
         DEC alpha, OSF1 3.x     (compile with: -Dalpha -DOSF1V3)
         DEC alpha, OSF1 3.x     (compile with: -Dalpha -DOSF1V4)

       --                                              -------------- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "error.h"

extern char   messageStr[];

#if !defined(SETNOTRAPS)

#if defined(mips)

/* --- Traps floating point exceptions for mips/IRIX5.

  CPU: mips
   OS: IRIX5

  see: ``man handle_sigfpes''

       Requires linking with -lfpe
       --                                              -------------- */

#include <limits.h>
#include <sigfpe.h>


void SetFPEtraps(void)
{
  sigfpe_[_UNDERFL].repls = _ZERO;

  sigfpe_[_OVERFL].abort  = sigfpe_[_OVERFL].trace  = 1;
  sigfpe_[_DIVZERO].abort = sigfpe_[_DIVZERO].trace = 1;
  sigfpe_[_INVALID].abort = sigfpe_[_INVALID].trace = 1;

  handle_sigfpes(_ON, _EN_OVERFL | _EN_DIVZERO | _EN_INVALID,
                  0, _ABORT_ON_ERROR, 0);
  Error(MESSAGE, "SetFPEtraps", 
	"\n-Setting FPE traps for mips (IRIX 5.x)\n");
}

/* ------- end ------------------------------------------------------ */

#elif defined(sparc)

#include <signal.h>

#if defined(SunOS4)

/* --- Traps floating point exceptions for sparc/SunOS4.

  CPU: sparc
   OS: SunOS4

  see: ``man ieee_handler''

       --                                              -------------- */

void  Trapped_FPE_Exception(int sig, int code, struct sigcontext *scp, 
		             char *address);
int   ieee_handler(const char *action, const char  *exception,
                   sigfpe_handler_type hdl);

void SetFPEtraps(void)
{
  const char routineName[] = "SetFPEtraps";

  if (ieee_handler("set", "common",
       (sigfpe_handler_type) Trapped_FPE_Exception) != 0)
    Error(MESSAGE, routineName, "IEEE trapping not supported here\n");
  else
    Error(MESSAGE, routineName,
	  "\n-Setting FPE traps for sparc (SunOS 4.x)\n");
}

void Trapped_FPE_Exception(int sig, int code, struct sigcontext *scp, 
		  char *address)
{
  const char routineName[] = "Trapped_FPE_Exception";
  char *type = "unknown";
 
  switch(code) {
 
  case FPE_INTDIV_TRAP:   type = "integer divide by zero          "; break;
  case FPE_INTOVF_TRAP:   type = "integer overflow                "; break;
 
  case FPE_FLTDIV_TRAP:   type = "floating point division by zero "; break;
  case FPE_FLTUND_TRAP:   type = "floating point underflow        "; break;
  case FPE_FLTOVF_TRAP:   type = "floating point overflow         "; break;
  case FPE_FLTINEX_TRAP:  type = "floating point inexact          "; break;
  }
 
  sprintf(messageStr, "  ---- trapped IEEE FPE: %s ----\n"
	  "  ---- signal: %d, code: 0x%x, addr: 0x%x, aborting ----\n",
	  type, sig, code, address);
  Error(MESSAGE, routineName, messageStr);
  abort();
}
/* ------- end ------------------------------------------------------ */


#else

/* --- Traps floating point exceptions for sparc/SunOS5.

  CPU: sparc
   OS: SunOS5

  see: ``man ieee_handler''

       Requires linking with -lsunmath
       --                                              -------------- */
 
#include <stdlib.h>
#include <sunmath.h>
#include <siginfo.h>
#include <ucontext.h>

void Trapped_FPE_Exception(int sig, siginfo_t *sip, ucontext_t *uap);

void SetFPEtraps(void)
{
  const char routineName[] = "SetFPEtraps";

  if (ieee_handler("set", "common", 
                   (sigfpe_handler_type) Trapped_FPE_Exception) != 0)
    Error(MESSAGE, routineName, "IEEE trapping not supported here\n");
  else
    Error(MESSAGE, routineName,
	  "\n-Setting FPE traps for sparc (SunOS 5.x)\n");
}

void Trapped_FPE_Exception(int sig, siginfo_t *sip, ucontext_t *uap)
{
  const char routineName[] = "Trapped_FPE_Exception";

  char *type = "unknown";
 
  switch(sip->si_code) {
 
  case FPE_INTDIV:   type = "integer divide by zero          "; break;
  case FPE_INTOVF:   type = "integer overflow                "; break;
  case FPE_FLTDIV:   type = "floating point division by zero "; break;
  case FPE_FLTUND:   type = "floating point underflow        "; break;
  case FPE_FLTOVF:   type = "floating point overflow         "; break;
  case FPE_FLTRES:   type = "floating point inexact          "; break;
  case FPE_FLTINV:   type = "invalid floating point operation"; break;
  case FPE_FLTSUB:   type = "subscript out of range          "; break;
  }
 
  sprintf(messageStr, "  ---- trapped IEEE FPE: %s ----\n"
	  "  ---- signal: %d, code: 0x%x, aborting ----\n",
	  type, sig, sip->si_code);
  Error(MESSAGE, routineName, messageStr);
  abort();
}
#endif

/* ------- end ------------------------------------------------------ */


#elif defined(alpha)

/* --- Traps floating point exceptions for DEC alpha/OSF1V[34]

  CPU: alpha
   OS: OSF1V[34]

  see: ``man ieee''

       --                                              -------------- */

#include <signal.h>
#include <machine/fpu.h>

void Trapped_FPE_Exception(int dummy);

void SetFPEtraps()
{
  const char routineName[] = "SetFPEtraps";

  unsigned long fp_control = 0;
  struct sigaction action, o_action;

  fp_control = (IEEE_TRAP_ENABLE_INV |
		IEEE_TRAP_ENABLE_DZE |
		IEEE_TRAP_ENABLE_OVF |
		IEEE_TRAP_ENABLE_UNF);
  ieee_set_fp_control(fp_control);
  action.sa_handler = Trapped_FPE_Exception;

  if (sigaction(SIGFPE, &action, &o_action) != 0)
    Error(MESSAGE, routineName, "\nIEEE trapping not supported here\n");
  else
    Error(MESSAGE, routineName,
	  "\nSetting FPE traps for DEC alpha (OSF1 3.x)\n");
}

void Trapped_FPE_Exception(int dummy)
{
  const char routineName[] = "Trapped_FPE_Exception";

  unsigned long fp_control = ieee_get_fp_control();

  sprintf(messageStr, "  ---- fp_control = 0x%lx\n", fp_control);
  Error(MESSAGE, messageStr);

  if (fp_control & IEEE_STATUS_INV)
    Error(MESSAGE, routineName,
	  "  ---- trapped IEEE FPE: Invalid operation\n");
  if (fp_control & IEEE_STATUS_DZE)
    Error(MESSAGE, routineName,
	  "  ---- trapped IEEE FPE: floating point division by zero\n");
  if (fp_control & IEEE_STATUS_INE)
    Error(MESSAGE, routineName,
	  "  ---- trapped IEEE FPE: floating point inexact\n");
  if (fp_control & IEEE_STATUS_UNF)
    Error(MESSAGE, routineName,
	  "  ---- trapped IEEE FPE: floating point underflow\n");
  if (fp_control & IEEE_STATUS_OVF)
    Error(MESSAGE, routineName, 
	  "  ---- trapped IEEE FPE: floating point overflow\n");

  abort();
}
/* ------- end ------------------------------------------------------ */

#else  /* if defined(mips, sparc, alpha) */

/* ------- unknown -------------------------------------------------- */

void SetFPEtraps(void)
{
  Error(MESSAGE, "SetFPEtraps",
	"\nUnsupported CPU and/or OS: cannot set FPE traps explicitly\n");
}

#endif

#else

void SetFPEtraps(void)
{
  /* --- Explicitly do not set traps --                -------------- */

  Error(MESSAGE, "SetFPEtraps", 
	"\nFPE traps have not been set explicitly\n");
}

#endif /* #if !defined(SETNOTRAPS) */
/* ------- end ------------------------------------------------------ */
