#ifndef LOCARNA_OPTIONS_HH
#define LOCARNA_OPTIONS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/*------------------------------------------------------------

  Copyright (C) 1999 by Sebastian Will.

  All Rights Reserved.
  
  ------------------------------------------------------------*/

/************************************************************
 *
 * options -- an interface for getopt_long
 *
 ************************************************************/

/************************************************************
 * 
 * interface to getopt
 *  
 *  - easier to use
 *  - may be restricted to certain cases
 *  - hopefully useable for many cases
 *
 ************************************************************/

#include <getopt.h>

#include <string>

namespace LocARNA {

#ifndef FALSE
#  define FALSE 0
#endif

#ifndef TRUE
#  define TRUE 1
#endif

    /* argument types */
#define O_NO_ARG 0
#define O_ARG_STRING 1
#define O_ARG_INT 2
#define O_ARG_FLOAT 3
#define O_ARG_DOUBLE 4
#define O_ARG_BOOL 5
#define O_SECTION -1
#define O_SECTION_HIDE -2

#define O_NODEFAULT std::string("__")

    /**
       @brief Definition structure of an option
    */
    typedef
    struct {
	std::string longname; //!<  long option name 
	char shortname; //!<  short option char 
	bool *flag; //!<  pointer to flag that indicates if option given 
	int arg_type; //!<  type of argument 
	void *argument; //!<  pointer to variable that should hold argument, 0 indicates no arg 
	std::string deflt; //!<  pointer to default argument, if arg optional. otherwise 0 
	std::string argname; //!< optional name for an argument (shown in usage string)
	std::string description; //!< optional description (shown in help) 
    }
	option_def;

    /* longname==0 and shortname==0 and arg_type<=O_SECTION is not allowed for regular options definition */

    /*
      Example for option_def array

      bool opt_help;
      .
      .
      .
      int optVal_size;
      int default_size=1000;

      struct option_def my_options[] = {
      {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,0,"This help"},
      {"num",'n',&opt_num,O_ARG_INT,&optVal_num,O_NODEFAULT,0,"Some arbitrary number"},
      {"output",'o',0,O_ARG_STRING,&outputfile,O_NODEFAULT,
      "output-file", "File for output"}, // mandatory
      {"size",'s',0,O_ARG_INT,&optVal_size,"100","size","Size of problem"},
      {0,0,0,O_ARG_STRING,&inputfile,O_NODEFAULT,"input-file","File for input"},
      {0,0,0,0,0,0,0,0}
      }; 

      the last entry shows how to define non-option command line
      arguments. If there is more than one such definition, the order of
      those argument definitions is important. This is different to the
      option argument definitions.

    */

    /* ***********************************************************/
    /* error message */
    extern char *O_error_msg;

    /* ***********************************************************/
    /* prototypes */

    /** process options */
    bool process_options(int argc, char *argv[], option_def options[]);

    /* print a usage string */
    void print_usage(char *progname, option_def options[]);

    /* print a longer help */
    void print_help(char *progname, option_def options[]);

    /** 
     * Print all options and their settings to standard out
     * 
     * @param options options array 
     */
    void print_options(option_def options[]);

} // end namespace LocARNA

#endif
