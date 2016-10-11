#ifndef LOCARNA_STOPWATCH_HH
#define LOCARNA_STOPWATCH_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <tr1/unordered_map>
#include <iosfwd>
#include <string>


namespace LocARNA {    
    /**
     * @brief control a set of named stop watch like timers
     */
    class StopWatch {
    private:
	struct timer_t {
	    bool running; //!<whether the timer is running
	    double last_start; //!< last start time
	    double total; //!< total accumulated time
	    size_t cycles; //!<number of start/stop cycles
	    
	    timer_t(): running(false), last_start(0.0), total(0.0), cycles(0) {}
	};

	//! type of map to store named timers
	typedef std::tr1::unordered_map<std::string,timer_t> map_t;
	
	map_t timers;
	
	bool print_on_exit;

    public:
	
	/** 
	 * @brief Constructor
	 * 
	 * @param print_on_exit whether to automatically print times on exit  
	 */
	StopWatch(bool print_on_exit=false);
	

	/** 
	 * @brief Destructor 
	 */
	~StopWatch();
	
	/** 
	 * Control automatic printing of times at exit
	 * 
	 * @param print_on_exit whether to print on exit
	 */
	void
	set_print_on_exit(bool print_on_exit);
	
	/** 
	 * @brief start a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return success
	 */
	bool
	start(const std::string &name);
	
	/** 
	 * @brief stop a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return success
	 */
	bool
	stop(const std::string &name);

	/** 
	 * @brief test whether named timer is running
	 * 
	 * @param name timer name
	 * 
	 * @return running?
	 */
	bool
	is_running(const std::string &name) const;

	/** 
	 * @brief current total time of a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return time (if running add time since start)
	 */
	double
	current_total(const std::string &name) const;
	
	/** 
	 * @brief current start/stop cycles of a named timer
	 * 
	 * @param name timer name
	 * 
	 * @return cycles (including started cycle if running)
	 */
	size_t current_cycles(const std::string &name) const;
	
	/** 
	 * @brief print information for one timer
	 * 
	 * @param out output stream
	 * @param name 
	 *
	 * @note determine current running time for running timers
	 *
	 * @return output stream
	 */
	std::ostream &
	print_info(std::ostream &out,const std::string &name) const;

	/** 
	 * @brief print information for all timers
	 * 
	 * @param out output stream
	 * 
	 * @return output stream
	 */
	std::ostream &
	print_info(std::ostream &out) const;

	
    private:
	double current_time () const;
    };
}
#endif
