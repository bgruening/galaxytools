use threads;
use threads::shared;
use Thread::Semaphore;



## execute sub for all given argument lists in in parallel
##
# @param $sub_ref       --- reference to the sub-routine
# @param $argument_lists_ref --- reference to the list of argument-lists references
# @param $thread_num     --- number of parallel threads
#
# idea: use $thread_num worker threads. each thread gets a new job as
# long as there is one and processes the job. job reservation needs to
# be synchronized (by lock on shared var $job_count)
#
sub foreach_par {
    my ($sub_ref,$argument_lists_ref,$thread_num) = @_;
    
    my @argument_lists = @{ $argument_lists_ref };
    
    my $job_num=$#argument_lists+1; ## total number of jobs

    my $job_count : shared = 0; ## count already started jobs

    my @my_threads;
    
    for (my $i=1; $i<=$thread_num; $i++) {
	$my_threads[$i]=
	    threads->create(
 		sub {
		    while(1) {
			# get the job $job_count
			my $my_job;
			{
			    lock($job_count);
			    $my_job=$job_count;
			    $job_count++;
			}
			
			## if job number too large then terminate thread
			if ($my_job >= $job_num) {return;}
			
			## start job
			$sub_ref->(@{ $argument_lists[$my_job] });
		    }
		}
	    );
    }

    for (my $i=1; $i<=$thread_num; $i++) {
 	$my_threads[$i]->join();
    }
}


## generate a thread-unique filename
sub threadsafe_name($) {
    my ($name) = @_;
    if (threads->tid()==0) { return $name; }
    else {return $name.(threads->tid());}
}


return 1;
