#ifdef _MSC_VER
#define DLLEXPORT(func) DEC$ ATTRIBUTES DLLEXPORT :: func
#define DLLEXPORT_ALIAS(func,alias) DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"alias" :: func
#else
#define DLLEXPORT(func)
#define DLLEXPORT_ALIAS(func,alias)  
#endif