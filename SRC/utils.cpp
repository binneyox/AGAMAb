#include "utils.h"
#include "math_core.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include <signal.h>
// stack trace presumably only works with GCC under unix-compatible platforms
#if(defined(__GNUC__) && !defined(__MINGW32__))
#define HAVE_STACKTRACE
#include <cxxabi.h>
#include <execinfo.h>
#endif

namespace utils {

EXP bool fileExists(const std::string& fileName)
{
    std::ifstream infile(fileName.c_str());
    return infile.good();
}

/* -------- error reporting routines ------- */

namespace{  // internal

/// remove function signature from GCC __PRETTY_FUNCTION__
inline std::string undecorateFunction(std::string origin)
{
#ifdef __GNUC__
    // parse the full function signature returned by __PRETTY_FUNCTION__
    std::string::size_type ind=origin.find('(');
    if(ind!=std::string::npos)
        origin.erase(ind);
    ind=origin.rfind(' ');
    if(ind!=std::string::npos)
        origin.erase(0, ind+1);
#endif
    return origin;
}

#ifdef HAVE_STACKTRACE
/// check if a character is a part of a valid C++ identifier name
inline bool isAlphanumeric(const char c)
{
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_';
}

/// attempt to demangle any possible C++ identifier names that are encountered in the input string
std::string demangleName(char* const input)
{
    std::string result;
    char *ptr = input, *name = NULL;
    while(1) {                       // scan the entire input line until '\0' is found.
        if(isAlphanumeric(*ptr)) {   // current character is a part of an identifier name:
            if(!name)                // if we are not already in the regime of recording a name,
                name = ptr;          // turn this regime on and remember the start of the name.
        } else {                     // otherwise the character is not part of a name:
            if(name) {               // if we have been recording a name, finalize this name.
                char oldsym = *ptr;  // remember the character at the current position
                *ptr = '\0';         // temporarily construct a null-terminated string with the name
                int status;          // try to demangle an identifier name
                char* demangledName = abi::__cxa_demangle(name, NULL, NULL, &status);
                if(status == 0) {    // success: append the demangled name to the results
                    result.append(demangledName);
                } else {             // demangling failed: append the name itself
                    result.append(name);
                }
                free(demangledName); // free the temporary string returned by __cxa_demangle()
                *ptr = oldsym;       // put the current character back into place
            }
            name = NULL;             // turn off the regime of recording a name
            if(*ptr != '\0')         // append the current character (not part of any name) to results
                result.append(1, *ptr);
        }
        if(*ptr == '\0')             // reached the end of the input string
            break;
        else
            ++ptr;                   // pass on to the next character of the input string
    }
    return result;
}
#endif

/// a prefix character prepended to the output message
inline char verbosityText(VerbosityLevel level)
{
    switch(level) {
        case VL_MESSAGE: return ':';
        case VL_WARNING: return '!';
        case VL_DEBUG:   return '-';
        case VL_VERBOSE: return '=';
        default: return '?';
    }
}

/// file to store the messages sent to msg() routine;
/// if not open, they are printed to stderr
static std::ofstream logfile;

/// read the environment variables controlling the verbosity level and log file redirection
VerbosityLevel initVerbosityLevel()
{
//	const char* env = std::getenv("LOGFILE");
	char *env;
	size_t len;
	_dupenv_s(&env,&len,"LOGFILE");
    if(len>1) {
        logfile.open(env);
    }
//    env = std::getenv("LOGLEVEL");
    _dupenv_s(&env,&len,"LOGLEVEL");
    if(env && env[0] >= '0' && env[0] <= '3')
        return static_cast<VerbosityLevel>(env[0]-'0');
    return VL_MESSAGE;  // default
}

/// default routine that dumps text messages to stderr
void defaultmsg(VerbosityLevel level, const char* origin, const std::string &message)
{
	if(level > verbosityLevel)
		return;
	std::string msg = verbosityText(level) + 
			  (origin==NULL || *origin=='\0' ? "" : "{"+undecorateFunction(origin)+"} ") + message + '\n';
	if(logfile.is_open())
		logfile << msg << std::flush;
	else
		std::cerr << msg << std::flush;
}

}  // namespace

/// global pointer to the routine that displays or logs information messages
//EXP MsgType* msg = &defaultmsg;

EXP void msg(VerbosityLevel level, const char* origin, const std::string &message)
{
	if(level > verbosityLevel)
		return;
	std::string msg = verbosityText(level) + 
			  (origin==NULL || *origin=='\0' ? "" : "{"+undecorateFunction(origin)+"} ") + message + '\n';
	if(logfile.is_open())
		logfile << msg << std::flush;
	else
		std::cerr << msg << std::flush;
}


/// global variable controlling the verbosity of printout
EXP VerbosityLevel verbosityLevel = initVerbosityLevel();
EXP VerbosityLevel getVerbosityLevel(void){return verbosityLevel;}

EXP std::string stacktrace()
{
#ifdef HAVE_STACKTRACE
    const int MAXFRAMES = 256;
    void* tmp[MAXFRAMES];     // temporary storage for stack frames
    int numframes = backtrace(tmp, MAXFRAMES);
    // convert this information into strings with function signatures and offsets
    char** lines = backtrace_symbols(tmp, numframes);  // this array should be freed later
    // represent these strings in a more user-friendly way (demangle identifier names)
    std::string result = "Stack trace ("+toString(numframes)+
        " functions, starting from the most recent one):\n";
    for(int i=1; i<numframes; i++)
        result += demangleName(lines[i]) + '\n';
    free(lines);
    return result;
#else
    return "Stack trace not available\n";
#endif
}

/* ----------- mechanism for capturing the Control-Break signal ----------------- */
namespace {
/// flag that is set to one if the keyboard interrupt has occurred
volatile bool ctrlBreakTriggered;
/// signal handler installed during lengthy computations that triggers the flag
void customCtrlBreakHandler(int) { ctrlBreakTriggered = true; }
/// a stack of previous signal handlers restored after the computation is finished
std::vector<void(*)(int)> prevCtrlBreakHandler;
}

EXP CtrlBreakHandler::CtrlBreakHandler()
{
    // this class could be instantiated multiple times in nested routines,
    // but the break flag is cleared only for the outermost one.
    if(prevCtrlBreakHandler.empty())
        ctrlBreakTriggered = false;
    // store the previous signal handler on the stack, and set the new one
    prevCtrlBreakHandler.push_back(signal(SIGINT, customCtrlBreakHandler));
}

EXP CtrlBreakHandler::~CtrlBreakHandler()
{
    // restore the previous handler once the instance of the class is destroyed
    assert(!prevCtrlBreakHandler.empty());       // it must have been set in the constructor
    signal(SIGINT, prevCtrlBreakHandler.back()); // restore the previous handler
    prevCtrlBreakHandler.pop_back();             // and eliminate it from the stack
}

EXP bool CtrlBreakHandler::triggered() { return ctrlBreakTriggered; }

/* ----------- string/number conversion and parsing routines ----------------- */

EXP int toInt(const char* val) {
    double result = toDouble(val);  // to allow input values such as "1e5"
    if(!isFinite(result))
        throw std::invalid_argument("Parse error: \"" + std::string(val) +
            "\" does not contain a valid integer number");
    if(result >= (double)std::numeric_limits<int>::max())
        return std::numeric_limits<int>::max();
    if(result <= (double)std::numeric_limits<int>::min())
        return std::numeric_limits<int>::min();
    return (int)trunc(result);
}

EXP float toFloat(const char* val) {
    if(*val == '\0')
        return 0;
    char* end;
    float result = strtof(val, &end);
    if(end == val)  // could not parse a single character
        throw std::invalid_argument("Parse error: \"" + std::string(val) +
            "\" does not contain a valid float number");
    return result;
}

EXP double toDouble(const char* val) {
    if(*val == '\0')
        return 0;
    char* end;
    double result = strtod(val, &end);
    if(end == val)  // could not parse a single character
        throw std::invalid_argument("Parse error: \"" + std::string(val) +
            "\" does not contain a valid double number");
    return result;
}

EXP bool toBool(const char* val) {
    while(1)   // scan the string until a valid character is encountered or an exception is thrown
    switch(*val) {
        case '\0':  // empty string counts as false
        case '0':
        case 'n':
        case 'N':
        case 'f':
        case 'F': return false;
        case '1':
        case 'y':
        case 'Y':
        case 't':
        case 'T': return true;
        case ' ':
        case '\t':  ++val; break;  // skip whitespace and continue from the next character
        default:  throw std::invalid_argument("Parse error: \"" + std::string(val) +
            "\" does not contain a valid boolean value");
    }
}

EXP std::string toString(double val, unsigned int width) {
    char buf[100];
    snprintf(buf, 100, "%-.*g", width, val);
    return std::string(buf);
}

EXP std::string toString(float val, unsigned int width) {
    char buf[100];
    snprintf(buf, 100, "%-.*g", width, val);
    return std::string(buf);
}

EXP std::string toString(int val) {
    char buf[100];
    snprintf(buf, 100, "%i", val);
    return std::string(buf);
}

EXP std::string toString(unsigned int val) {
    char buf[100];
    snprintf(buf, 100, "%u", val);
    return std::string(buf);
}

EXP std::string toString(long val) {
    char buf[100];
    snprintf(buf, 100, "%li", val);
    return std::string(buf);
}

EXP std::string toString(unsigned long val) {
    char buf[100];
    snprintf(buf, 100, "%lu", val);
    return std::string(buf);
}

EXP std::string toString(long long val) {
    char buf[100];
    snprintf(buf, 100, "%lli", val);
    return std::string(buf);
}

EXP std::string toString(unsigned long long val) {
    char buf[100];
    snprintf(buf, 100, "%llu", val);
    return std::string(buf);
}

EXP std::string toString(const void* val) {
    char buf[100];
    snprintf(buf, 100, "%p", val);
    return std::string(buf);
}

//  Pretty-print - convert float (and integer) numbers to string of fixed width.
//  Employ sophisticated techniques to fit the number into a string of exactly the given length.
EXP std::string pp(double num, unsigned int uwidth)
{
    const int MAXWIDTH = 31;  // rather arbitrary restriction, but doubles never are that long
    std::string result;
    int width = std::min<int>(uwidth, MAXWIDTH);
    if(width<1)
        return result;
    bool sign = std::signbit(num);
    if(num==0) {
        result.resize(width, '0');
        if(width>1 && sign) result[0] = '-';
        if(width>1+sign) result[1+sign] = '.';
        return result;
    }
    if(num!=num || num==INFINITY || num==-INFINITY) {
        result = num==INFINITY ? "+INF" : num==-INFINITY ? "-INF" : "NAN";
        result.resize(width, ' ');
        return result;
    }
    // separate out sign, and reduce the available width
    if(sign) {
        result = "-";
        width--;
    }
    if(width==0) {
        return result;
    }
    num = fabs(num);

    // decimal exponent
    int expon = (int)floor(log10(num));
    char buf[MAXWIDTH+1];

    // now we have three possibilities:
    // 1) exponential notation:  2.34e5 or 6e-7 - at least one digit before the 'e' symbol
    // 2) fixed-point notation:  234321.2 or 0.0000006321
    // 3) "#" if none of the above fits into the available space

    // try to print the number x>=1 in fixed-point format
    if(expon >= 0 && expon <= width-1) {
        int len = snprintf(buf, MAXWIDTH, "%-#.*f", std::max<int>(width-2-expon, 0), num);
        if(len == width+1 && buf[width] == '.') {
            // exceeds the required width, but the final decimal point may be removed
            buf[width] = '\0';
            len--;
        }
        if(len == width) {
            result += buf;
            return result;
        }
        // otherwise we may have encountered the situation of rounding up
        // a number 9.x to 10. and exceeding the width
    }

    // expected length of the exponent part of the string if we choose exponential notation
    // (including the 'e' symbol, and removing the sign of exponent if it is positive)
    int  len_exp = expon<=-100 ? 5 : expon<=-10 ? 4 : expon<0 ? 3 : expon<10 ? 2 : expon<100 ? 3 : 4;

    // expected # of significant digits in mantissa in the exp format (including the leading digit)
    int  dig_exp = std::max<int>(width-len_exp-1, 1);  // (one position is occupied by decimal point)

    // expected number of significant digits in fixed-point format if |x|<1
    // (first two symbols are "0.", and then possibly a few more zeros before the number begins)
    int  dig_lt1 = width-1+expon;  // e.g. if expon=-2 and width=6, we may have 0.0234 - 3 digit accuracy

    // try to print the number x<1 in fixed-point format if the expected precision is no less than
    // in the exponential format, or in the special case of a number 0.5<=x<1 rounded up to 1.
    if(expon < 0 && (dig_lt1 >= dig_exp || (num>=0.5 && width<=2))) {
        int len = snprintf(buf, MAXWIDTH, "%-#.*f", std::max<int>(width-2, 0), num);
        if(len > width) {
            buf[width] = '\0';
            len--;
        }
        result += buf;
        return result;
    }

    // try to print the number in exponential format
    if(width >= len_exp+1) {
        // strip out exponent, so that the number is within [1:10)
        num /= math::pow(10., expon);
        if(num<1 || !isFinite(num))  // it might be <1 due to roundoff error
            num = 1.;
        int len = snprintf(buf, MAXWIDTH, "%-#.*f", std::max<int>(width-2-len_exp, 0), num);
        if(len >= 2 && buf[0] == '1' && buf[1] == '0') {
            // a number 9.x is rounded up to 10, so we should replace it with 1. and increase the exponent
            expon++;
            buf[1] = '.';
            if(len>2)
                buf[2] = '0';
        }
        char buf_exp[8] = {'e'};
        len_exp = snprintf(buf_exp+1, 6, "%-i", expon)+1;
        if(len_exp < width) {
            if(len > width-len_exp)
                buf[width-len_exp] = '\0';
            result += buf;
            result += buf_exp;
            return result;
        }
    }

    // can't use any of these - display a jail symbol
    result.resize(sign+width, '#');
    return result;
}

EXP std::vector<std::string> splitString(const std::string& src, const std::string& delim)
{
    std::vector<std::string> result;
    std::string str(src);
    std::string::size_type indx=str.find_first_not_of(delim);
    if(indx==std::string::npos) {
        result.push_back("");   // ensure that result contains at least one element
        return result;
    }
    if(indx>0)  // remove unnecessary delimiters at the beginning
        str=str.erase(0, indx);
    while(!str.empty()) {
        indx=str.find_first_of(delim);
        if(indx==std::string::npos)
            indx=str.size();
        result.push_back(str.substr(0, indx));
        str=str.erase(0, indx);
        indx=str.find_first_not_of(delim);
        if(indx==std::string::npos)
            break;
        str=str.erase(0, indx);
    }
    return result;
}

EXP bool endsWithStr(const std::string& str, const std::string& end)
{
    return end.size()<=str.size() && str.find(end, str.size()-end.size())!=str.npos;
}

EXP bool stringsEqual(const std::string& str1, const std::string& str2)
{
    std::string::size_type len=str1.size();
    if(len!=str2.size())
        return false;
    for(std::string::size_type i=0; i<len; i++)
        if(tolower((unsigned char)str1[i]) != tolower((unsigned char)str2[i]))
            return false;
    return true;
}

EXP bool stringsEqual(const std::string& str1, const char* str2)
{
    if(str2==NULL)
        return false;
    for(std::string::size_type i=0; i<str1.size(); i++)
        if(str2[i]==0 || tolower((unsigned char)str1[i]) != tolower((unsigned char)str2[i]))
            return false;
    return str2[str1.size()]==0;  // ensure that the 2nd string length is the same as the 1st
}

#define SIGN(A,B) (std::signbit(B)?(-A):(A))
EXP std::string jjbcmp(double x){
	const int i1 = 90, i2 = 8100, i3 = 729000;
	const double at2t23 = 8388608., big=1.e28, small=1.e-28, aln2=log(2.);
	char m[6]={80,0,0,0,0,'\0'};
	int lnx,rem; double absx=fabs(x);
	if (x!=0){
		if (absx>big) x= SIGN(big, x);
		if (absx<small) x = SIGN(small, x);
		if (absx>=1.)
			lnx = (int)(log(absx)/aln2) + 129;
		else
			lnx = (int)(log(absx)/aln2) + 128;
		int j = lnx/i1;
		m[3] = lnx - j*i1;
		lnx = lnx - 129;
		if (lnx>0)
			rem = (int)((absx/pow(2.,(double)(lnx)) - 1.)*at2t23);
		else
			rem = (int)((absx*pow(2.,(double)(fabs(lnx))) - 1.)*at2t23);
		m[4] = rem/i3;
		rem = rem - m[4]*i3;
		m[2] = rem/i2;
		rem = rem - m[2]*i2;
		m[1] = rem/i1;
		m[0] = rem - m[1]*i1;
		m[4] = m[4] + j*12;
		if (x<0.) m[4] = m[4] + 40;
	}
	for (int i=0;i<5;i++){
		m[i] = m[i] + 33;
		if (m[i]>=94) m[i] = m[i] + 1;
		if (m[i]>=123) m[i] = m[i] + 1;
	}
	std::string st(m);
	return st;
}

EXP bool jjbget(FILE *tmpf,float *xg,int npt){
	float expon, at2t23= 8388608.;
	long int i1=90;
	unsigned char m[5];
	bool ok=true;
	for(int k=0; k<npt; k++){
		m[0]=getc(tmpf); if(m[0]=='\n') m[0]=getc(tmpf);
		if(feof(tmpf)){
			ok=false; break;
		}

		for(int j=1;j<5;j++) m[j]=getc(tmpf);
		for(int i = 0; i<5; i++){
			if (m[i]>=124) m[i] = m[i] - 1;
			if (m[i]>=95) m[i] = m[i] - 1;
			m[i] = m[i] - 33;
		}
		if (m[4]==80) xg[k] = 0;
		else{
			if (m[4]>=40){
				m[4] = m[4] - 40;
				xg[k] = -1.;
			}else xg[k] = 1;
			int j = m[4]/12;
			m[4] = m[4] - j*12;
			expon=(float)(j*i1+m[3]-129);
			xg[k]=xg[k]*(float)((((m[4]*i1+m[2])*i1+m[1])*i1
					     +m[0])/at2t23+ 1.);
			if (expon>0.)
				xg[k]=xg[k]*pow(2.,expon);
			else
				xg[k]=xg[k]/pow(2.,-expon);
		}
	}
	return ok;
}

}  // namespace
