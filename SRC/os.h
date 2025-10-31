//To get visual studio to export symbols from .dll library

//If have visual studio export symbols from .dll library
#ifdef _MSC_VER
#define EXP __declspec(dllexport)
#else
#define _dupenv_s(env,len,filename){\
        (*env)=getenv(filename);\
        (*len)=sizeof(*env);\
}
#define fscanf_s fscanf
#define EXP
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen(filename,mode))==NULL
#endif
