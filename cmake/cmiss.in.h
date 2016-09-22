/*
* Header file configured by CMake to enable builds for various platforms
*/
#cmakedefine HAVE_SA_NODEFER 

#ifndef HAVE_SA_NODEFER
#define	SA_NODEFER 0x40000000u
#endif

#cmakedefine HAVE_SIGACTION_STRUCT
#ifndef HAVE_SIGACTION_STRUCT
struct sigaction {
	int dummy;
	void* sa_handler;
	unsigned long sa_flags;
	/*__sigrestore_t sa_restorer;*/
	int sa_mask;
	/*sigset_t sa_mask;*/		/* mask last for extensibility */
};
#define sigaction(a,b,c) a
#define sigemptyset(a) a
#endif